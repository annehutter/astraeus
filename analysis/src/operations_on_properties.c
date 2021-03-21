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
#include "derive_properties.h"
#include "get_grid_properties.h"
#include "map_gal_to_grid.h"
#include "get_further_derived_properties.h"
#include "operations_on_properties.h"

int check_operation_present(char *propertyName)
{
  int result = 0;
  if(strchr(propertyName, '*') != NULL || strchr(propertyName, '%') != NULL)
    result = 1;
  return result;
}

void get_delimiters_in_propertyName(char *propertyName, int *numDelim, char **listDelim)
{
  int tmpNumDelim = 10;
  char *tmpListDelim = allocate_array_char(tmpNumDelim, "tmpListDelim");
  char *tmpPropertyName = propertyName;
  char *tmp1 = NULL, *tmp2 = NULL;
  
  int i = 0;
  while(strchr(tmpPropertyName, '*') != NULL || strchr(tmpPropertyName, '%') != NULL)
  {
    if(strchr(tmpPropertyName, '*') != NULL && strchr(tmpPropertyName, '%') != NULL)
    {
      tmp1 = strchr(tmpPropertyName, '*') + 1;
      tmp2 = strchr(tmpPropertyName, '%') + 1;
      if(strlen(tmp1) > strlen(tmp2))
      {
        tmpPropertyName = tmp1;
        tmpListDelim[i] = '*';
      }
      else
      {
        tmpPropertyName = tmp2;
        tmpListDelim[i] = '%';
      }
    }
    else if(strchr(tmpPropertyName, '*') != NULL)
    {
      tmpPropertyName = strchr(tmpPropertyName, '*') + 1;
      tmpListDelim[i] = '*';
    }
    else if(strchr(tmpPropertyName, '%') != NULL)
    {
      tmpPropertyName = strchr(tmpPropertyName, '%') + 1;
      tmpListDelim[i] = '%';
    }
    i++;
    if(i >= tmpNumDelim)
    {
      tmpListDelim = realloc(tmpListDelim, (int)(tmpNumDelim * 1.2));
    }
  }
  *numDelim = i;
  
  tmpListDelim = realloc(tmpListDelim, (*numDelim)*sizeof(char));
  if(tmpListDelim == NULL)
  {
    fprintf(stderr, "Could not realloc tmpListDelim\n");
    exit(EXIT_FAILURE);
  }
  *listDelim = tmpListDelim;
}

void extract_properties(char *propertyName, int *numProperties, char ***listPropertyNames)
{
  char copyPropertyName[strlen(propertyName)+1];
  strcpy(copyPropertyName, propertyName);
  char *tmpPropertyName = propertyName;
  int numDelim = 0;
  char *listDelim = NULL;
  get_delimiters_in_propertyName(copyPropertyName, &numDelim, &listDelim);

  int maxOffset = strlen(propertyName);
  int offset = 0;
  char **tmpListPropertyNames = malloc(sizeof(char*)*(numDelim+1));
  char *token = NULL;
  
  int i = 0;
  char thisString[2];
  thisString[1] = '\0';
  
  while(offset < maxOffset && i<numDelim)
  {
    tmpPropertyName = &(propertyName[offset]);
    thisString[0] = listDelim[i];
    token = strtok(tmpPropertyName, thisString);
    offset += strlen(token) + 1;
    tmpListPropertyNames[i] = token;
    i++;
  }
  tmpListPropertyNames[i] = &(propertyName[offset]);

  *numProperties = i+1;
  *listPropertyNames = tmpListPropertyNames;
  
  free(listDelim);
}

void multiply_properties(int num, double *property1, double* property2)
{
  for(int i=0; i<num; i++)
    if(property1[i] != -1. && property2[i] != -1.)
      property1[i] = property1[i]*property2[i];
}

void divide_properties(int num, double *property1, double* property2)
{
  for(int i=0; i<num; i++)
    if(property1[i] != -1. && property2[i] != -1.)
      property1[i] = property1[i]/property2[i];
}

int check_MvirThreshold_present(char *propertyName)
{
  int result = 0;
  if(strstr(propertyName, "_MVIRCUT") != NULL)
    result = 1;
  
  return result;
}

char *get_propertyName(char *propertyName)
{
  int length = strlen(propertyName);
  char *newPropertyName = NULL;
    
  if(strstr(propertyName, "_MVIRCUT") != NULL)
  {
    newPropertyName = malloc(sizeof(char)*(length - 8 + 1));
    newPropertyName = memcpy(newPropertyName, propertyName, sizeof(char)*(length-8));
    memcpy(newPropertyName+length-8, propertyName+length, sizeof(char));
  }
  else
  {
    newPropertyName = malloc(sizeof(char)*(length + 1));
    newPropertyName = memcpy(newPropertyName, propertyName, sizeof(char)*(length+1));
  }
      
  return newPropertyName;
}

int get_num_reducedProperty(int numGal, int numSnaps, double *Mvir, double MvirThreshold)
{
  int counter = 0;
  
  for(int gal=0; gal<numGal; gal++)
  {
    if(Mvir[gal * numSnaps + numSnaps - 1] >= MvirThreshold)
      counter++;
  }
  
  return counter;
}

double *reduceProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times, double *thisProperty, int *numGalSelectedProperty)
{
  int numGal = 0;
  double *Mvir = getThisProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "Mvir", currSnap, times, &numGal);
  double MvirThreshold = simParam->MvirThreshold;
    
  int numReducedGal = get_num_reducedProperty(numGal, 1, Mvir, MvirThreshold);
  double *reducedProperty = allocate_array_double(numReducedGal, "reducedProperty");
  
  int counter = 0;
  for(int gal=0; gal<numGal; gal++)
  {
    if(Mvir[gal] >= MvirThreshold)
    {
      reducedProperty[counter] = thisProperty[gal];
      counter++;
    }
  }
  
  *numGalSelectedProperty = numReducedGal;
  
  free(Mvir);
  
  return reducedProperty;
}

double *reducePropertyHistory(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times, double *thisPropertyHistory, int *numGalSelectedProperty)
{
  int numGal = 0;
  int numSnaps = currSnap + 1;
  double *Mvir = getThisProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "Mvir", currSnap, times, &numGal);
  double MvirThreshold = simParam->MvirThreshold;
  
  int numReducedGal = get_num_reducedProperty(numGal, 1, Mvir, MvirThreshold);
  double *reducedPropertyHistory = allocate_array_double(numReducedGal * numSnaps, "reducedProperty");
  
  int counter = 0;
  for(int gal=0; gal<numGal; gal++)
  {
    if(Mvir[gal] >= MvirThreshold)
    {
      for(int snap=0; snap<numSnaps; snap++)
        reducedPropertyHistory[counter * numSnaps + snap] = thisPropertyHistory[gal * numSnaps + snap];
      counter++;
    }
  }
  
  *numGalSelectedProperty = numReducedGal;
  
  free(Mvir);
  
  return reducedPropertyHistory;
}

double *selectPropertyWithMapping(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, double propertyLowLimit, double propertyUpLimit, int currSnap, float *times)
{
  double *newProperty = NULL;
  
  /* check whether property includes operation */
  if(check_operation_present(propertyName) == 0)
  {
    if(check_grid_property_present(propertyName) > 0)
    {
      newProperty = getGridProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times);
    }
    else if(check_map_gal_to_grid_property_present(propertyName) > 0)
    {
      newProperty = get_map_gal_to_grid_property(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, propertyLowLimit, propertyUpLimit, currSnap, times);
    }
    else if(check_derived_property_present(propertyName) > 0)
    {
      newProperty = get_derived_property(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times);
    }
    else
    {
      newProperty = getProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times);
    }
  }
  
  return newProperty;
}

double *selectProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  double *newProperty = NULL;
  
  /* check whether property includes operation */
  if(check_operation_present(propertyName) == 0)
  {
    if(check_grid_property_present(propertyName) > 0)
    {
      newProperty = getGridProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times);
    }
    else if(check_derived_property_present(propertyName) > 0)
    {
      newProperty = get_derived_property(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times);
    }
    else
    {
      newProperty = getProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times);
    }
  }
  
  return newProperty;
}

double *getThisPropertyWithMapping(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, double propertyLowLimit, double propertyUpLimit, int currSnap, float *times, int *numGalProperty)
{
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
  double *newProperty = NULL;
  double *newOtherProperty = NULL;
  char *newPropertyName = get_propertyName(propertyName);
  char copyPropertyName[strlen(newPropertyName)+1];
  strcpy(copyPropertyName, newPropertyName);

  /* check whether property includes operation */
  if(check_operation_present(newPropertyName) == 0)
  {
    newProperty = selectPropertyWithMapping(simParam, treeList, numTrees, listEquals, index, sizeListEquals, newPropertyName, propertyLowLimit, propertyUpLimit, currSnap, times);
  }
  else
  {    
    int numDelim = 0;
    char *listDelim = NULL;
    get_delimiters_in_propertyName(copyPropertyName, &numDelim, &listDelim);
    
    int numProperties = 0;
    char **listPropertyNames;
    extract_properties(copyPropertyName, &numProperties, &listPropertyNames);
    
    newProperty = selectProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[0], currSnap, times);

    char thisString[2];
    thisString[1] = '\0';
    for(int i=0; i<numDelim; i++)
    {
      newOtherProperty = selectProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[i+1], currSnap, times);
      
      thisString[0] = listDelim[i];
      if(strcmp(thisString, "*") == 0)
        multiply_properties(numGal, newProperty, newOtherProperty);
      else if(strcmp(thisString, "%") == 0)
        divide_properties(numGal, newProperty, newOtherProperty);
    }
    
    free(listDelim);
    free(listPropertyNames);
    free(newOtherProperty);
  }
  
  *numGalProperty = numGal;
    
  /* check whether property includes lower cut in Mvir */
  if(check_MvirThreshold_present(propertyName) == 1)
  {
    int numGalReducedProperty = 0;
    double *newReducedProperty = reduceProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, currSnap, times, newProperty, &numGalReducedProperty);
    
    *numGalProperty = numGalReducedProperty;
    if(newProperty != NULL) free(newProperty);
    newProperty = newReducedProperty;
    newReducedProperty = NULL;
  }
  
  if(newPropertyName != NULL) free(newPropertyName);
  
  return newProperty;
}

double *getThisProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times, int *numGalProperty)
{
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
  double *newProperty = NULL;
  double *newOtherProperty = NULL;
  char *newPropertyName = get_propertyName(propertyName);
  char copyPropertyName[strlen(newPropertyName)+1];
  strcpy(copyPropertyName, newPropertyName);
      
  /* check whether property includes operation */
  if(check_operation_present(newPropertyName) == 0)
  {
    newProperty = selectProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, newPropertyName, currSnap, times);
  }
  else
  {    
    int numDelim = 0;
    char *listDelim = NULL;
    get_delimiters_in_propertyName(copyPropertyName, &numDelim, &listDelim);
    
    int numProperties = 0;
    char **listPropertyNames;
    extract_properties(copyPropertyName, &numProperties, &listPropertyNames);
    
    newProperty = selectProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[0], currSnap, times);

    char thisString[2];
    thisString[1] = '\0';
    for(int i=0; i<numDelim; i++)
    {
      newOtherProperty = selectProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[i+1], currSnap, times);
      
      thisString[0] = listDelim[i];
      if(strcmp(thisString, "*") == 0)
        multiply_properties(numGal, newProperty, newOtherProperty);
      else if(strcmp(thisString, "%") == 0)
        divide_properties(numGal, newProperty, newOtherProperty);
    }
    
    free(listDelim);
    free(listPropertyNames);
    free(newOtherProperty);
  }
  
  *numGalProperty = numGal;
    
  /* check whether property includes lower cut in Mvir */
  if(check_MvirThreshold_present(propertyName) == 1)
  {
    int numGalReducedProperty = 0;
    double *newReducedProperty = reduceProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, currSnap, times, newProperty, &numGalReducedProperty);
    
    *numGalProperty = numGalReducedProperty;
    if(newProperty != NULL) free(newProperty);
    newProperty = newReducedProperty;
    newReducedProperty = NULL;
  }
  
  if(newPropertyName != NULL) free(newPropertyName);
  
  return newProperty;
}

double *getThisPropertyHistory(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times, int *numGalProperty)
{
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
  double *newProperty = NULL;
  double *newOtherProperty = NULL;
  char *newPropertyName = get_propertyName(propertyName);
  char copyPropertyName[strlen(newPropertyName)+1];
  strcpy(copyPropertyName, newPropertyName);

  /* check whether property includes operation */
  if(check_operation_present(newPropertyName) == 0)
  {
    newProperty = getPropertyHistory(simParam, treeList, numTrees, listEquals, index, sizeListEquals, newPropertyName, currSnap, times);
  }
  else
  {
    int numSnaps = currSnap + 1;
    
    int numDelim = 0;
    char *listDelim = NULL;
    get_delimiters_in_propertyName(copyPropertyName, &numDelim, &listDelim);
    
    int numProperties = 0;
    char **listPropertyNames;
    extract_properties(copyPropertyName, &numProperties, &listPropertyNames);
    
    newProperty = getPropertyHistory(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[0], currSnap, times);

    char thisString[2];
    thisString[1] = '\0';
    for(int i=0; i<numDelim; i++)
    {
      newOtherProperty = getPropertyHistory(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[i+1], currSnap, times);
      
      thisString[0] = listDelim[i];
      if(strcmp(thisString, "*") == 0)
        multiply_properties(numGal * numSnaps, newProperty, newOtherProperty);
      else if(strcmp(thisString, "%") == 0)
        divide_properties(numGal * numSnaps, newProperty, newOtherProperty);
    }
    
    free(listDelim);
    free(listPropertyNames);
    free(newOtherProperty);
  }
  
  *numGalProperty = numGal;
    
  /* check whether property includes lower cut in Mvir */
  if(check_MvirThreshold_present(propertyName) == 1)
  {
    int numGalReducedProperty = 0;
    double *newReducedProperty = reducePropertyHistory(simParam, treeList, numTrees, listEquals, index, sizeListEquals, currSnap, times, newProperty, &numGalReducedProperty);
        
    *numGalProperty = numGalReducedProperty;
    if(newProperty != NULL) free(newProperty);
    newProperty = newReducedProperty;
    newReducedProperty = NULL;
  }
  
  if(newPropertyName != NULL) free(newPropertyName);
  
  return newProperty;
}