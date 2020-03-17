#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "dconfObj.h"
#include "outgal.h"
#include "zas_lists.h"

#include "build_index_tree_walking.h"
#include "derive_properties.h"
#include "operations_on_properties.h"
// #include "statistics.h"
// #include "selection.h"
#include "statistics_evolution.h"
#include "selection_evolution.h"

void get_1D_histogram_evolution(dconfObj_t simParam, outgtree_t **theseTrees, int numTrees, int **index, int ***listEquals, int *sizeListEquals, int **listMerged, int listID, int thisRank)
{    
  int binsInLog = simParam->binsInLog_1D_evolution[listID];
  int binsPerMag = simParam->binsPerMag_1D_evolution[listID];
  float binsMinValue = simParam->binsMinValue_1D_evolution[listID];
  float binsMaxValue = simParam->binsMaxValue_1D_evolution[listID];
    
  int numBins1 = get_numBins_for_fixed_bins(binsMinValue, binsMaxValue, binsInLog, binsPerMag);
  int minIndex1 = get_minIndex_for_fixed_bins(binsMinValue, binsMaxValue, binsInLog, binsPerMag);
  float *values1 = get_values_for_fixed_bins(binsMinValue, binsMaxValue, binsInLog, binsPerMag);  

  /* allocate histogram */
  int numSnaps = simParam->numSnaps;
  float *histogram = allocate_array_float(numBins1 * numSnaps, "histogram");
  int *histogramNum = allocate_array_int(numBins1 * numSnaps, "histogramNum");

  int currSnap = 0, prevSnap = -1;
  for(int snap=0; snap<numSnaps; snap++)
  { 
    if(thisRank == 0) printf("snap = %d\n", snap);
    
    currSnap = snap;
    get_listEquals(theseTrees, numTrees, index, listEquals, sizeListEquals, listMerged, prevSnap, currSnap);
    prevSnap = currSnap;
    
    /* reading data */
    float *binProperty1 = NULL;
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, simParam->binProperty_1D_evolution[listID], currSnap, simParam->times);
    
    float *property = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, simParam->property_1D_evolution[listID], currSnap, simParam->times);
    
    /* fill histogram with data */
    if(get_numGal_at_snap(sizeListEquals, numTrees) > 0)
      calc_1D_histogram_evolution(get_numGal_at_snap(sizeListEquals, numTrees), property, currSnap, binProperty1, binsInLog, binsPerMag, minIndex1, 0, numBins1, numSnaps, &histogram, &histogramNum);

    free(property);
    free(binProperty1);
  }
  
  float volumeInMpc = pow(simParam->boxsize/simParam->hubble_h, 3);
  float factor = 1./volumeInMpc;
  
  for(int valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(int snap=0; snap<numSnaps; snap++)
    {
      histogram[valueInt1 * numSnaps + snap] = factor * histogram[valueInt1 * numSnaps + snap];
    }
  }
  
  /* communicate histograms */
#ifdef MPI
  float *recvHistogram = allocate_array_float(numBins1 * numSnaps, "recvHistogram");
  int *recvNumHistogram = allocate_array_int(numBins1 * numSnaps, "recvNumHistogram");
  
  MPI_Allreduce(histogram, recvHistogram, numBins1*numSnaps, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1*numSnaps, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(int valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(int snap=0; snap<numSnaps; snap++)
    {
      histogram[valueInt1 * numSnaps + snap] = recvHistogram[valueInt1 * numSnaps + snap];
      histogramNum[valueInt1 * numSnaps + snap] = recvNumHistogram[valueInt1 * numSnaps + snap];
    }
  }
  
  free(recvHistogram);
  free(recvNumHistogram);
#endif
  
  /* save histogram */
  char *filename = create_filename_1D_histogram_evolution(simParam->outputDir, simParam->property_1D_evolution[listID], simParam->binProperty_1D_evolution[listID]);
  write_1D_histogram_evolution(numBins1, values1, numSnaps, simParam->redshifts, histogram, histogramNum, thisRank, filename);
  free(filename);

  /* deallocate histogram */
  free(values1);
  free(histogram);
  free(histogramNum);
}

char *create_filename_1D_histogram_evolution(char *directory, char *property, char *binProperty)
{
  char *filename = concat_strings(6, directory, "/1Dhistogram_evolution_", property, "_", binProperty, ".dat");
  
  return filename;
}