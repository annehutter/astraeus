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
#include "statistics.h"

#include "selection.h"
#include "operations_on_properties.h"

void get_2D_histogram_history(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, char *filename, int thisRank)
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
  
  calc_2D_histogram_history(numGal, property, currSnap+1, simParam->redshifts, currSnap, binProperty1, binProperty2, binsInLog1, binsInLog2, binsPerMag1, binsPerMag2, containsHistories1, containsHistories2, thisRank, filename);
  
  free(property);
  free(binProperty1);
  free(binProperty2);
}

void get_1D_histogram_history(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, int binsInLog1, int binsPerMag1, int containsHistories1, char *filename, int thisRank)
{
  int numGal1 = 0;
  double *binProperty1 = NULL;
  
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  else
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  
  int numGal = 0;
  double *property = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, simParam->times, &numGal);
  
  assert(numGal1 == numGal);
  
  calc_1D_histogram_history(numGal, property, currSnap+1, simParam->redshifts, currSnap, binProperty1, binsInLog1, binsPerMag1, containsHistories1, thisRank, filename);
  
  free(property);
  free(binProperty1);
}

void get_3D_histogram_value(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, char *binProperty2Name, char *binProperty3Name, int binsInLog1, int binsInLog2, int binsInLog3, int binsPerMag1, int binsPerMag2, int binsPerMag3, int containsHistories1, int containsHistories2, int containsHistories3, double binProperty1LowLimit, double binProperty1UpLimit, char *filename, int thisRank)
{
  int numGal1 = 0, numGal2 = 0, numGal3 = 0;
  double *binProperty1 = NULL, *binProperty2 = NULL, *binProperty3 = NULL;
    
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  else
    binProperty1 = getThisPropertyWithMapping(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, binProperty1LowLimit, binProperty1UpLimit, currSnap, simParam->times, &numGal1);
  
  if(containsHistories2 == 1)
    binProperty2 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times, &numGal2);
  else
    binProperty2 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty2Name, currSnap, simParam->times, &numGal2);
  
  if(containsHistories3 == 1)
    binProperty3 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty3Name, currSnap, simParam->times, &numGal3);
  else
    binProperty3 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty3Name, currSnap, simParam->times, &numGal3);
  
  int numGal = 0;
  double *property = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, simParam->times, &numGal);
  
  assert(numGal1 == numGal);
  assert(numGal2 == numGal);
  assert(numGal3 == numGal);
    
  calc_3D_histogram_value(numGal, property, currSnap, binProperty1, binProperty2, binProperty3, binsInLog1, binsInLog2, binsInLog3, binsPerMag1, binsPerMag2, binsPerMag3, containsHistories1, containsHistories2, containsHistories3, thisRank, filename);
  
  free(property);
  free(binProperty1);
  free(binProperty2);
  free(binProperty3);
}

void get_2D_histogram_value(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, double binProperty1LowLimit, double binProperty1UpLimit, char *filename, int thisRank)
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
  
  calc_2D_histogram_value(numGal, property, currSnap, binProperty1, binProperty2, binsInLog1, binsInLog2, binsPerMag1, binsPerMag2, containsHistories1, containsHistories2, thisRank, filename);
  
  free(property);
  free(binProperty1);
  free(binProperty2);
}

void get_1D_histogram_value(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, int binsInLog1, int binsPerMag1, int containsHistories1, double binProperty1LowLimit, double binProperty1UpLimit, char *filename, int thisRank)
{
  int numGal1 = 0;
  double *binProperty1 = NULL;
  
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  else
    binProperty1 = getThisPropertyWithMapping(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, binProperty1LowLimit, binProperty1UpLimit, currSnap, simParam->times, &numGal1);
  
  int numGal = 0;
  double *property = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, simParam->times, &numGal);
  
  assert(numGal1 == numGal);
  
  calc_1D_histogram_value(numGal, property, currSnap, binProperty1, binsInLog1, binsPerMag1, containsHistories1, thisRank, filename);
  
  free(property);
  free(binProperty1);
}

void get_2D_histogram(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, char *filename, int thisRank)
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
      
  assert(numGal1 == numGal2);
  
  float volumeInMpc = pow(simParam->boxsize/simParam->hubble_h, 3);
  
  calc_2D_histogram(numGal1, binProperty1, binProperty2, currSnap, binsInLog1, binsInLog2, binsPerMag1, binsPerMag2, containsHistories1, containsHistories2, volumeInMpc, thisRank, filename);
  
  free(binProperty1);
  free(binProperty2);
}

void get_1D_histogram(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *binProperty1Name, int binsInLog1, int binsPerMag1, int containsHistories1, int cumulative, char *filename, int thisRank)
{
  int numGal1  = 0;
  double *binProperty1 = NULL;
  
  if(containsHistories1 == 1)
    binProperty1 = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  else
    binProperty1 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, binProperty1Name, currSnap, simParam->times, &numGal1);
  
  float volumeInMpc = pow(simParam->boxsize/simParam->hubble_h, 3);
  
  calc_1D_histogram(numGal1, binProperty1, currSnap, binsInLog1, binsPerMag1, containsHistories1, volumeInMpc, cumulative, thisRank, filename);
  
  free(binProperty1);
}

char *create_filename_2D_histogram(char *directory, char *property, char *binProperty1, char *binProperty2, float redshift, float smoothingScale, double MvirThreshold)
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
    filename = concat_strings(11, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, smoothingScale_name, redshift_name, ".dat");
  }
  else if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL)
  {
    filename = concat_strings(10, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, redshift_name, ".dat");
  }
  else if(strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL)
  {
    filename = concat_strings(10, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, smoothingScale_name, redshift_name, ".dat");
  }
  else
  {
    filename = concat_strings(9, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, redshift_name, ".dat");
  }
  
  free(propertyName);
  free(binProperty1Name);
  free(binProperty2Name);
  
  return filename;
}

char *create_filename_1D_histogram(char *directory, char *property, char *binProperty, float redshift, float smoothingScale, double MvirThreshold)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);
  
  char smoothingScale_name[12];
  if(strstr(binProperty, "smooth") != NULL || strstr(binProperty, "SMOOTH") != NULL)
  {
    sprintf(smoothingScale_name, "%3.2f", smoothingScale);
  }

  char *propertyName = NULL;
  char *binPropertyName = NULL;
  if(strstr(property, "_MVIRCUT") != NULL)
  {
    char Mvir_name[12];
    sprintf(Mvir_name, "%.2e", MvirThreshold);
    
    propertyName = concat_strings(2, property, Mvir_name);
    binPropertyName = concat_strings(2, binProperty, Mvir_name);
  }
  else
  {
    propertyName = concat_strings(1, property);
    binPropertyName = concat_strings(1, binProperty);
  }
  
  char *filename = NULL;
  
  if(strstr(binProperty, "smooth") != NULL || strstr(binProperty, "SMOOTH") != NULL)
  {
    filename = concat_strings(8, directory, "/1Dhistogram_", propertyName, "_", binPropertyName, smoothingScale_name, redshift_name, ".dat");
  }
  else
  {
    filename = concat_strings(7, directory, "/1Dhistogram_", propertyName, "_", binPropertyName, redshift_name, ".dat");
  }
  
  free(propertyName);
  free(binPropertyName);
  
  return filename;
}

char *create_filename_3D_histogram_value(char *directory, char *property, char *binProperty1, char *binProperty2, char *binProperty3, float redshift, float smoothingScale, double MvirThreshold)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);
  
  char smoothingScale_name[12];
  if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL || strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL || strstr(binProperty3, "smooth") != NULL || strstr(binProperty3, "SMOOTH") != NULL)
  {
    sprintf(smoothingScale_name, "%3.2f", smoothingScale);
  }

  char *propertyName = NULL;
  char *binProperty1Name = NULL;
  char *binProperty2Name = NULL;
  char *binProperty3Name = NULL;
  if(strstr(property, "_MVIRCUT") != NULL)
  {
    char Mvir_name[12];
    sprintf(Mvir_name, "%.2e", MvirThreshold);
    
    propertyName = concat_strings(2, property, Mvir_name);
    binProperty1Name = concat_strings(2, binProperty1, Mvir_name);
    binProperty2Name = concat_strings(2, binProperty2, Mvir_name);
    binProperty3Name = concat_strings(2, binProperty3, Mvir_name);
  }
  else
  {
    propertyName = concat_strings(1, property);
    binProperty1Name = concat_strings(1, binProperty1);
    binProperty2Name = concat_strings(1, binProperty2);
    binProperty3Name = concat_strings(1, binProperty3);
  }
  
  char *filename = NULL;
  
  if((strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL) && (strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL) && (strstr(binProperty3, "smooth") != NULL || strstr(binProperty3, "SMOOTH") != NULL))
  {
    filename = concat_strings(14, directory, "/3Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, smoothingScale_name, "-", binProperty3, smoothingScale_name, redshift_name, "_value.dat");
  }
  else if((strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL) && (strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL))
  {
    filename = concat_strings(13, directory, "/3Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, smoothingScale_name, "-", binProperty3, redshift_name, "_value.dat");
  }
  else if((strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL) && (strstr(binProperty3, "smooth") != NULL || strstr(binProperty3, "SMOOTH") != NULL))
  {
    filename = concat_strings(13, directory, "/3Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, "-", binProperty3, smoothingScale_name, redshift_name, "_value.dat");
  }
  else if((strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL) && (strstr(binProperty3, "smooth") != NULL || strstr(binProperty3, "SMOOTH") != NULL))
  {
    filename = concat_strings(13, directory, "/3Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, smoothingScale_name, "-", binProperty3, smoothingScale_name, redshift_name, "_value.dat");
  }
  else if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL)
  {
    filename = concat_strings(12, directory, "/3Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, "-", binProperty3, redshift_name, "_value.dat");
  }
  else if(strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL)
  {
    filename = concat_strings(12, directory, "/3Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, smoothingScale_name, "-", binProperty3, redshift_name, "_value.dat");
  }
  else if(strstr(binProperty3, "smooth") != NULL || strstr(binProperty3, "SMOOTH") != NULL)
  {
    filename = concat_strings(12, directory, "/3Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, "-", binProperty3, smoothingScale_name, redshift_name, "_value.dat");
  }
  else
  {
    filename = concat_strings(11, directory, "/3Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, "-", binProperty3, redshift_name, "_value.dat");
  }
  
  free(propertyName);
  free(binProperty1Name);
  free(binProperty2Name);
  free(binProperty3Name);
  
  return filename;
}

char *create_filename_2D_histogram_value(char *directory, char *property, char *binProperty1, char *binProperty2, float redshift, float smoothingScale, double MvirThreshold)
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
    filename = concat_strings(11, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, smoothingScale_name, redshift_name, "_value.dat");
  }
  else if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL)
  {
    filename = concat_strings(10, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, redshift_name, "_value.dat");
  }
  else if(strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL)
  {
    filename = concat_strings(10, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, smoothingScale_name, redshift_name, "_value.dat");
  }
  else
  {
    filename = concat_strings(9, directory, "/2Dhistogram_", propertyName, "_", binProperty1Name, "-", binProperty2Name, redshift_name, "_value.dat");
  }
  
  free(propertyName);
  free(binProperty1Name);
  free(binProperty2Name);
  
  return filename;
}

char *create_filename_1D_histogram_value(char *directory, char *property, char *binProperty, float redshift, float smoothingScale, double MvirThreshold)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);

  char smoothingScale_name[12];
  if(strstr(binProperty, "smooth") != NULL || strstr(binProperty, "SMOOTH") != NULL)
  {
    sprintf(smoothingScale_name, "%3.2f", smoothingScale);
  }

  char *propertyName = NULL;
  char *binPropertyName = NULL;
  if(strstr(property, "_MVIRCUT") != NULL)
  {
    char Mvir_name[12];
    sprintf(Mvir_name, "%.2e", MvirThreshold);
    
    propertyName = concat_strings(2, property, Mvir_name);
    binPropertyName = concat_strings(2, binProperty, Mvir_name);
  }
  else
  {
    propertyName = concat_strings(1, property);
    binPropertyName = concat_strings(1, binProperty);
  }
  
  char *filename = NULL;
  
  if(strstr(binProperty, "smooth") != NULL || strstr(binProperty, "SMOOTH") != NULL)
  {
    filename = concat_strings(8, directory, "/1Dhistogram_", propertyName, "_", binPropertyName, smoothingScale_name, redshift_name, "_value.dat");
  }
  else
  {
    filename = concat_strings(7, directory, "/1Dhistogram_", propertyName, "_", binPropertyName, redshift_name, "_value.dat");
  }
  
  free(propertyName);
  free(binPropertyName);
  
  return filename;
}

char *create_filename_2D_histogram_numDens(char *directory, char *binProperty1, char *binProperty2, float redshift, float smoothingScale, double MvirThreshold)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);

  char smoothingScale_name[12];
  if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL || strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL)
  {
    sprintf(smoothingScale_name, "%3.2f", smoothingScale);
  }

  char *binProperty1Name = NULL;
  char *binProperty2Name = NULL;
  if(strstr(binProperty1, "_MVIRCUT") != NULL)
  {
    char Mvir_name[12];
    sprintf(Mvir_name, "%.2e", MvirThreshold);
    
    binProperty1Name = concat_strings(2, binProperty1, Mvir_name);
    binProperty2Name = concat_strings(2, binProperty2, Mvir_name);
  }
  else
  {
    binProperty1Name = concat_strings(1, binProperty1);
    binProperty2Name = concat_strings(1, binProperty2);
  }
  
  char *filename = NULL;
  
  if((strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL) && (strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL))
  {
    filename = concat_strings(9, directory, "/2Dhistogram_numDens_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, smoothingScale_name, redshift_name, ".dat");
  }
  else if(strstr(binProperty1, "smooth") != NULL || strstr(binProperty1, "SMOOTH") != NULL)
  {
    filename = concat_strings(8, directory, "/2Dhistogram_numDens_", binProperty1Name, smoothingScale_name, "-", binProperty2Name, redshift_name, ".dat");
  }
  else if(strstr(binProperty2, "smooth") != NULL || strstr(binProperty2, "SMOOTH") != NULL)
  {
    filename = concat_strings(8, directory, "/2Dhistogram_numDens_", binProperty1Name, "-", binProperty2Name, smoothingScale_name, redshift_name, ".dat");
  }
  else
  {
    filename = concat_strings(7, directory, "/2Dhistogram_numDens_", binProperty1Name, "-", binProperty2Name, redshift_name, ".dat");
  }
  
  free(binProperty1Name);
  free(binProperty2Name);
  
  return filename;
}

char *create_filename_1D_histogram_numDens(char *directory, char *property, float redshift, float smoothingScale, double MvirThreshold)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);
  
  char smoothingScale_name[12];
  if(strstr(property, "smooth") != NULL || strstr(property, "SMOOTH") != NULL)
  {
    sprintf(smoothingScale_name, "%3.2f", smoothingScale);
  }

  char *propertyName = NULL;
  if(strstr(property, "_MVIRCUT") != NULL)
  {
    char Mvir_name[12];
    sprintf(Mvir_name, "%.2e", MvirThreshold);
    
    propertyName = concat_strings(2, property, Mvir_name);
  }
  else
  {
    propertyName = concat_strings(1, property);
  }
  
  char *filename = NULL;
  
  if(strstr(property, "smooth") != NULL || strstr(property, "SMOOTH") != NULL)
  {
    filename = concat_strings(6, directory, "/1Dhistogram_numDens_", propertyName, smoothingScale, redshift_name, ".dat");
  }
  else
  {
    filename = concat_strings(5, directory, "/1Dhistogram_numDens_", propertyName, redshift_name, ".dat");
  }
  
  free(propertyName);
  
  return filename;
}