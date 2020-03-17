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
#include "derive_properties.h"
#include "get_grid_properties.h"
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
  char copyPropertyName[strlen(propertyName)];
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

void multiply_properties(int num, float *property1, float* property2)
{
  for(int i=0; i<num; i++)
    if(property1[i] != -1. && property2[i] != -1.)
      property1[i] = property1[i]*property2[i];
}

void divide_properties(int num, float *property1, float* property2)
{
  for(int i=0; i<num; i++)
    if(property1[i] != -1. && property2[i] != -1.)
      property1[i] = property1[i]/property2[i];
}

float *getGridProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  char fileEnding[10];
  char *gridFilename = NULL;
  float *newProperty = NULL;
  int doublePrecision = 0;
  
  if(check_grid_property_present(propertyName) == 1)
  {
    sprintf(fileEnding, "_%03d", currSnap);
    gridFilename = concat(simParam->densFilename, fileEnding);
    doublePrecision = simParam->densInputInDoublePrecision;
  }
  else
  {
    sprintf(fileEnding, "_%02d", currSnap);
    gridFilename = concat(simParam->ionFilename, fileEnding);
    doublePrecision = simParam->ionInputInDoublePrecision;
  }
    
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
  float *posIndex = getProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "POS", currSnap, times);
  newProperty = get_grid_property(simParam->size, simParam->thisRank, gridFilename, numGal, posIndex, simParam->gridsize, doublePrecision, simParam->memoryIntensive);
  free(posIndex);
  free(gridFilename);
  
  return newProperty;
}

float *getThisProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  float *newProperty = NULL;
  float *newOtherProperty = NULL;
  char copyPropertyName[strlen(propertyName)];
  strcpy(copyPropertyName, propertyName);
  
  /* check whether property includes operation */
  if(check_operation_present(propertyName) == 0)
  {
    if(check_grid_property_present(propertyName) > 0)
    {
      newProperty = getGridProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times);
    }
    else
      newProperty = getProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times);
  }
  else
  {
    int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
    
    int numDelim = 0;
    char *listDelim = NULL;
    get_delimiters_in_propertyName(copyPropertyName, &numDelim, &listDelim);
    
    int numProperties = 0;
    char **listPropertyNames;
    extract_properties(copyPropertyName, &numProperties, &listPropertyNames);
    
    if(check_grid_property_present(propertyName) > 0)
    {
      newProperty = getGridProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[0], currSnap, times);
    }
    else
      newProperty = getProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[0], currSnap, times);

    char thisString[2];
    thisString[1] = '\0';
    for(int i=0; i<numDelim; i++)
    {
      if(check_grid_property_present(propertyName) > 0)
      {
        newOtherProperty = getGridProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[i+1], currSnap, times);
      }
      else
        newOtherProperty = getProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, listPropertyNames[i+1], currSnap, times);
      
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
  
  return newProperty;
}

float *getThisPropertyHistory(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  float *newProperty = NULL;
  float *newOtherProperty = NULL;
  char copyPropertyName[strlen(propertyName)];
  strcpy(copyPropertyName, propertyName);
  
  /* check whether property includes operation */
  if(check_operation_present(propertyName) == 0)
  {
    newProperty = getPropertyHistory(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times);
  }
  else
  {
    int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
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
  
  return newProperty;
}