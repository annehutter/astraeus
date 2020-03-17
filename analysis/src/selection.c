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
#include "read_trees.h"
#include "zas_lists.h"

#include "build_index_tree_walking.h"
#include "derive_properties.h"
#include "statistics.h"

#include "selection.h"
#include "operations_on_properties.h"

void get_2D_histogram_history(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, char *filename, int thisRank)
{
  float *binProperty1 = NULL, *binProperty2 = NULL;
  
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times);
  else
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times);
  
  if(containsHistories2 == 1)
    binProperty2 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times);
  else
    binProperty2 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times);
  
  float *property = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, simParam->times);
  
//   for(int snap=0; snap<=currSnap; snap++)
//     printf("rank %d: %d: %e\n", thisRank, snap, property[snap]);
  
  calc_2D_histogram_history(get_numGal_at_snap(sizeListEquals, numTrees), property, currSnap+1, simParam->redshifts, currSnap, binProperty1, binProperty2, binsInLog1, binsInLog2, binsPerMag1, binsPerMag2, containsHistories1, containsHistories2, thisRank, filename);
  
  free(property);
  free(binProperty1);
  free(binProperty2);
}

void get_1D_histogram_history(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, int binsInLog1, int binsPerMag1, int containsHistories1, char *filename, int thisRank)
{
  float *binProperty1 = NULL;
  
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times);
  else
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times);
  
  float *property = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, simParam->times);
  
  calc_1D_histogram_history(get_numGal_at_snap(sizeListEquals, numTrees), property, currSnap+1, simParam->redshifts, currSnap, binProperty1, binsInLog1, binsPerMag1, containsHistories1, thisRank, filename);
  
  free(property);
  free(binProperty1);
}

void get_2D_histogram(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, char *filename, int thisRank)
{
  float *binProperty1 = NULL, *binProperty2 = NULL;
  
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times);
  else
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times);
  
  if(containsHistories2 == 1)
    binProperty2 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times);
  else
    binProperty2 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times);
      
  float volumeInMpc = pow(simParam->boxsize/simParam->hubble_h, 3);
  
  calc_2D_histogram(get_numGal_at_snap(sizeListEquals, numTrees), binProperty1, binProperty2, currSnap, binsInLog1, binsInLog2, binsPerMag1, binsPerMag2, containsHistories1, containsHistories2, volumeInMpc, thisRank, filename);
  
  free(binProperty1);
  free(binProperty2);
}

void get_1D_histogram(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *binProperty1Name, int binsInLog1, int binsPerMag1, int containsHistories1, int cumulative, char *filename, int thisRank)
{
  float *binProperty1 = NULL;
  
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times);
  else
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times);
  
  float volumeInMpc = pow(simParam->boxsize/simParam->hubble_h, 3);
  
  calc_1D_histogram(get_numGal_at_snap(sizeListEquals, numTrees), binProperty1, currSnap, binsInLog1, binsPerMag1, containsHistories1, volumeInMpc, cumulative, thisRank, filename);
  
  free(binProperty1);
}

char *create_filename_2D_histogram(char *directory, char *property, char *binProperty1, char *binProperty2, float redshift)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);

  char *filename = concat_strings(9, directory, "/2Dhistogram_", property, "_", binProperty1, "-", binProperty2, redshift_name, ".dat");
  
  return filename;
}

char *create_filename_1D_histogram(char *directory, char *property, char *binProperty, float redshift)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);

  char *filename = concat_strings(7, directory, "/1Dhistogram_", property, "_", binProperty, redshift_name, ".dat");
  
  return filename;
}

char *create_filename_2D_histogram_numDens(char *directory, char *binProperty1, char *binProperty2, float redshift)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);

  char *filename = concat_strings(7, directory, "/2Dhistogram_numDens_", binProperty1, "-", binProperty2, redshift_name, ".dat");
  
  return filename;
}

char *create_filename_1D_histogram_numDens(char *directory, char *property, float redshift)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);

  char *filename = concat_strings(5, directory, "/1Dhistogram_numDens_", property, redshift_name, ".dat");
  
  return filename;
}