#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <assert.h>

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
#include "statistics_analytic.h"
#include "numdens_analytic.h"

#include "selection_analytic.h"
#include "operations_on_properties.h"

void get_2D_histogram_analytic(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, char *filename, char *hmfFilename, int thisRank)
{
  int numGal1 = 0, numGal2 = 0;
  double *binProperty1 = NULL, *binProperty2 = NULL, *halomassEndSnap = NULL, *numDensEndSnap = NULL;
  
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  else
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  
  if(containsHistories2 == 1)
    binProperty2 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times, &numGal2);
  else
    binProperty2 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times, &numGal2);
      
  assert(numGal1 == numGal2);
  
  halomassEndSnap = getProperty_endSnap(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Mvir", currSnap, simParam->times);
  
  numDensEndSnap = get_numDens(numGal1, halomassEndSnap, simParam->mergertreeBinwidth, hmfFilename, simParam->hubble_h);

  calc_2D_histogram_analytic(numGal1, binProperty1, binProperty2, numDensEndSnap, currSnap, binsInLog1, binsInLog2, binsPerMag1, binsPerMag2, containsHistories1, containsHistories2, thisRank, filename);
  
  free(binProperty1);
  free(binProperty2);
  free(halomassEndSnap);
  free(numDensEndSnap);
}

void get_1D_histogram_analytic(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *binProperty1Name, int binsInLog1, int binsPerMag1, int containsHistories1, int cumulative, char *filename, char *hmfFilename, int thisRank)
{
  int numGal1  = 0;
  double *binProperty1 = NULL, *halomassEndSnap = NULL, *numDensEndSnap = NULL;
  
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  else
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
 
  halomassEndSnap = getProperty_endSnap(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Mvir", currSnap, simParam->times);
  
  numDensEndSnap = get_numDens(numGal1, halomassEndSnap, simParam->mergertreeBinwidth, hmfFilename, simParam->hubble_h);
    
  calc_1D_histogram_analytic(numGal1, binProperty1, numDensEndSnap, currSnap, binsInLog1, binsPerMag1, containsHistories1, cumulative, thisRank, filename);
  
  free(binProperty1);
  free(halomassEndSnap);
  free(numDensEndSnap);
}