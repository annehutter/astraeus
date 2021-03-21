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
#include "operations_on_properties.h"
#include "statistics_median.h"

#include "selection_median.h"

void get_2D_histogram_history_median(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, char *filename, int thisRank, int size)
{
  int numGal1 = 0, numGal2 = 0;
  double *binProperty1 = NULL, *binProperty2 = NULL;
    
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  else
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  
  if(containsHistories2 == 1)
    binProperty2 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times, &numGal2);
  else
    binProperty2 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times, &numGal2);
  
  int numGal = 0;
  double *property = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, simParam->times, &numGal);
    
  assert(numGal1 == numGal);
  assert(numGal2 == numGal);
  
  calc_2D_histogram_history_median(numGal, property, currSnap+1, simParam->redshifts, currSnap, binProperty1, binProperty2, binsInLog1, binsInLog2, binsPerMag1, binsPerMag2, containsHistories1, containsHistories2, thisRank, size, filename);
  
  free(property);
  free(binProperty1);
  free(binProperty2);
}

void get_2D_histogram_median(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, double binProperty1LowLimit, double binProperty1UpLimit, char *filename, int thisRank, int size)
{
  int numGal1 = 0, numGal2 = 0;
  double *binProperty1 = NULL, *binProperty2 = NULL;
    
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  else
    binProperty1 = getThisPropertyWithMapping(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, binProperty1LowLimit, binProperty1UpLimit, currSnap, simParam->times, &numGal1);
  
  if(containsHistories2 == 1)
    binProperty2 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times, &numGal2);
  else
    binProperty2 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times, &numGal2);
  
  int numGal = 0;
  double *property = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, simParam->times, &numGal);
    
  assert(numGal1 == numGal);
  assert(numGal2 == numGal);
  
  calc_2D_histogram_median(numGal, property, currSnap, binProperty1, binProperty2, binsInLog1, binsInLog2, binsPerMag1, binsPerMag2, containsHistories1, containsHistories2, thisRank, size, filename);
  
  free(property);
  free(binProperty1);
  free(binProperty2);
}

char *create_filename_2D_histogram_history_median(char *directory, char *property, char *binProperty1, char *binProperty2, float redshift, float smoothingScale, double MvirThreshold)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);
  
  char *propertyName = NULL;
  char *binProperty1Name = NULL;
  char *binProperty2Name = NULL;
  if(strstr(property, "_MVIRCUT") != NULL)
  {
    char Mvir_name[12];
    sprintf(Mvir_name, "%.2e", MvirThreshold);
    
    propertyName = concat_strings(2, property, Mvir_name);
    binProperty1Name = concat_strings(2, binProperty1, Mvir_name);
    binProperty2Name = concat_strings(2, binProperty2, Mvir_name);
  }
  else
  {
    propertyName = concat_strings(1, property);
    binProperty1Name = concat_strings(1, binProperty1);
    binProperty2Name = concat_strings(1, binProperty2);
  }
  
  char smoothingScale_name[12];
  if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL || strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL)
  {
    sprintf(smoothingScale_name, "%3.2f", smoothingScale);
  }

  char *filename = NULL;
  
  if((strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL) && (strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL))
  {
    filename = concat_strings(11, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, smoothingScale_name, redshift_name, "_median.dat");
  }
  else if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL)
  {
    filename = concat_strings(10, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, redshift_name, "_median.dat");
  }
  else if(strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL)
  {
    filename = concat_strings(10, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, smoothingScale_name, redshift_name, "_median.dat");
  }
  else
  {
    filename = concat_strings(9, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, redshift_name, "_median.dat");
  }
  
  free(propertyName);
  free(binProperty1Name);
  free(binProperty2Name);
  
  return filename;
}

char *create_filename_2D_histogram_median(char *directory, char *property, char *binProperty1, char *binProperty2, float redshift, float smoothingScale, double MvirThreshold)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);
  
  char smoothingScale_name[12];
  if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL || strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL)
  {
    sprintf(smoothingScale_name, "%3.2f", smoothingScale);
  }

  char *propertyName = NULL;
  char *binProperty1Name = NULL;
  char *binProperty2Name = NULL;
  if(strstr(property, "_MVIRCUT") != NULL)
  {
    char Mvir_name[12];
    sprintf(Mvir_name, "%.2e", MvirThreshold);
    
    propertyName = concat_strings(2, property, Mvir_name);
    binProperty1Name = concat_strings(2, binProperty1, Mvir_name);
    binProperty2Name = concat_strings(2, binProperty2, Mvir_name);
  }
  else
  {
    propertyName = concat_strings(1, property);
    binProperty1Name = concat_strings(1, binProperty1);
    binProperty2Name = concat_strings(1, binProperty2);
  }
  
  char *filename = NULL;
  
  if((strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL) && (strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL))
  {
    filename = concat_strings(11, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, smoothingScale_name, redshift_name, "_median.dat");
  }
  else if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL)
  {
    filename = concat_strings(10, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, redshift_name, "_median.dat");
  }
  else if(strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL)
  {
    filename = concat_strings(10, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, smoothingScale_name, redshift_name, "_median.dat");
  }
  else
  {
    filename = concat_strings(9, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, redshift_name, "_median.dat");
  }
  
  free(propertyName);
  free(binProperty1Name);
  free(binProperty2Name);
  
  return filename;
}