#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"
#include "dconfObj.h"
#include "outgal.h"

#include "derive_properties.h"
#include "communicate_mpi.h"
#include "get_grid_properties.h"
#include "smoothing_grid.h"

#include "map_gal_to_grid.h"

int check_map_gal_to_grid_property_present(char *propertyName)
{
  int result = 0;
  
  /* check whether a map to gal grid property is being asked for */
  if(strstr(propertyName, "GRID") != NULL)
    result = 1;
  
  return result;
}

char *get_map_gal_to_grid_propertyName(char *propertyName)
{
  int length = strlen(propertyName);
  char *newPropertyName = NULL;
  
  if(strstr(propertyName, "SMOOTH") != NULL)
  {
    newPropertyName = malloc(sizeof(char)*(length - 10 + 1));
    newPropertyName = memcpy(newPropertyName, propertyName+4, sizeof(char)*(length-10));
    memcpy(newPropertyName+length-10, propertyName+length, sizeof(char));
  }
  else
  {
    newPropertyName = malloc(sizeof(char)*(length - 4 + 1));
    newPropertyName = memcpy(newPropertyName, propertyName+4, sizeof(char)*(length-4+1));
  }
      
  return newPropertyName;
}

double *get_mapped_gal_positions(int numGal, double *galProperty, double *posIndex, double lowLimit, double upLimit, int *selectedNumGal)
{  
  double *selectedPosIndex = allocate_array_double(numGal, "selectedPosIndex");
  int counter = 0;
  
  for(int gal=0; gal<numGal; gal++)
  {
    if(galProperty[gal] >= lowLimit && galProperty[gal] < upLimit)
    {
      selectedPosIndex[counter] = posIndex[gal];
      counter++;
    }
  }
  
  *selectedNumGal = counter;
  selectedPosIndex = realloc(selectedPosIndex, sizeof(double)*counter);
  
  return selectedPosIndex;
}

#ifdef MPI
float *create_grid_with_map_gal_mpi(int size, int thisRank, int numEntries, double *posIndex, int gridsize)
{
  /*--------------------------------*/
  /* prepare all necessary arrays */
  /*--------------------------------*/

  int *indexArray = NULL;
  
  /* create rankArray */
  int *rankArray = create_rankArray_3Dgrid(numEntries, posIndex, gridsize);

  /* create sortedArray = posIndex sorted by ranks */
  double *sortedArray = create_sortedArray_double(numEntries, posIndex, rankArray, &indexArray);
  
  /* create array with number of galaxies to send to ranks */
  int32_t *numOnRanks = create_numToSendToRanksArray(size, numEntries, rankArray);
  
  /* arrays for receiving */
  int32_t *recvNumOnRanks = allocate_array_int32_t(size, "recvNumOnRanks");
  int32_t numRecvEntries = 0;
  double *recvArray = NULL;
    
  /* grid parameters */
  float *thisGrid = NULL;
  int zNbins = 0;
  int offset = 0;
  int gridIndex = 0;
  
  /*--------------------------------*/
  /* sending forward */
  /*--------------------------------*/  
  
  send_recv_array_double(size, thisRank, numOnRanks, &recvNumOnRanks, sortedArray, &recvArray);  
  free(sortedArray);
  sortedArray = NULL;
  
  /*--------------------------------*/
  /* map galaxies to grid */
  /*--------------------------------*/  
  
  numRecvEntries = get_numEntries_numOnThisRank(size, recvNumOnRanks);
  
  /* create grid */
  get_grid_decomposition(gridsize, &zNbins, &offset, 0);
  thisGrid = allocate_array_float(gridsize * gridsize * zNbins, "thisGrid");
  
  /* mapping galaxies to grid */
  for(int gal=0; gal<numRecvEntries; gal++)
  {
    gridIndex = (int)(recvArray[gal]) - offset*gridsize*gridsize;
    if(gridIndex < 0 || gridIndex >= gridsize*gridsize*zNbins)
      printf("rank %d: gal = %d\tgridIndex = %d\t offset = %d\t zNbins = %d\tPos = %f %d\n", thisRank, gal, gridIndex, offset, zNbins, recvArray[gal], (int)(recvArray[gal]));
    thisGrid[gridIndex] = thisGrid[gridIndex] + 1;
  }

  /*--------------------------------*/
  /* deallocation */
  /*--------------------------------*/
  free(indexArray);
  free(rankArray);
  free(sortedArray);
  free(numOnRanks);
  free(recvNumOnRanks);
  free(recvArray);

  return thisGrid;
}
#endif

#ifdef MPI
float *create_grid_with_map_gal_mpi_memintensive(int numEntries, double *posIndex, int gridsize)
{
  float *localGrid = allocate_array_float(gridsize * gridsize * gridsize, "thisGrid");
  float *thisGrid = allocate_array_float(gridsize * gridsize * gridsize, "thisGrid");
  int gridIndex = 0;
  
  /* mapping galaxies to grid */
  for(int gal=0; gal<numEntries; gal++)
  {
    gridIndex = (int)(posIndex[gal]);
    localGrid[gridIndex] = localGrid[gridIndex] + 1;
  }
  
  MPI_Allreduce(localGrid, thisGrid, gridsize * gridsize *gridsize, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  
  free(localGrid);

  return thisGrid;
}
#endif

float *create_grid_with_map_gal_serial(int numEntries, double *posIndex, int gridsize)
{
  float *thisGrid = allocate_array_float(gridsize*gridsize*gridsize, "thisGrid");
  int gridIndex = 0;
  
  for(int gal=0; gal<numEntries; gal++)
  {
    gridIndex = (int)(posIndex[gal]);
    thisGrid[gridIndex] = thisGrid[gridIndex] + 1;
  }
  
  return thisGrid;
}

float *create_grid_with_map_gal(int size, int thisRank, int numEntries, double *posIndex, int gridsize, int memoryIntensive)
{
  float *thisGrid = NULL;
    
#ifdef MPI
  if(memoryIntensive == 0)
  {
    thisGrid = create_grid_with_map_gal_mpi(size, thisRank, numEntries, posIndex, gridsize);
  }
  else
  {
    thisGrid = create_grid_with_map_gal_mpi_memintensive(numEntries, posIndex, gridsize);
  }
#else
  thisGrid = create_grid_with_map_gal_serial(numEntries, posIndex, gridsize);
#endif

  return thisGrid;
}

double *get_map_gal_to_grid_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, double galPropertyLowLimit, double galPropertyUpLimit, int currSnap, float *times)
{
  double smoothingScale = (simParam->smoothingScale / simParam->boxsize) * simParam->gridsize;
  double *thisProperty = NULL;
  
  /* get which galaxy property to load */
  char *galPropertyName = get_map_gal_to_grid_propertyName(propertyName);
    
  /* load galaxy property and positions */
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees); 
  double *galProperty = getProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, galPropertyName, currSnap, times);
  double *posIndex = getProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "POS", currSnap, times);

  /* get positions of selected galaxies */
  int selectedNumGal = 0;
  double *selectedPosIndex = get_mapped_gal_positions(numGal, galProperty, posIndex, galPropertyLowLimit, galPropertyUpLimit, &selectedNumGal);
    
  free(galPropertyName);
  free(galProperty);
  
  if(selectedNumGal > 0)
  {
    /* create grid */
    float *thisGrid = create_grid_with_map_gal(simParam->size, simParam->thisRank, selectedNumGal, selectedPosIndex, simParam->gridsize, simParam->memoryIntensive);
    free(selectedPosIndex);

    if(strstr(propertyName, "SMOOTH") != NULL)
    {
      smooth_grid_property(&thisGrid, simParam->gridsize, simParam->memoryIntensive, smoothingScale);
    }
    
    /* get grid value for galaxies */
#ifdef MPI
    if(simParam->memoryIntensive == 0)
    {
      thisProperty = communicate_array_double(simParam->size, simParam->thisRank, numGal, posIndex, simParam->gridsize, posIndex, thisGrid);
    }
    else
    {
      thisProperty = allocate_array_double(numGal, "thisProperty");
      get_grid_values(numGal, &thisProperty, posIndex, simParam->gridsize, thisGrid, 1);
    }
#else
    thisProperty = allocate_array_double(numGal, "thisProperty");
    get_grid_values(numGal, &thisProperty, posIndex, simParam->gridsize, thisGrid, 1);
#endif
    
    free(thisGrid);
  }
  
  free(posIndex);
  
  return thisProperty;
}
