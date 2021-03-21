#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

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
#include "smoothing_grid.h"
#include "get_grid_properties.h"

/*----------------------------------------------------------------*/
/* DO DOMAIN DECOMPOSITION FOR GRIDS */
/*----------------------------------------------------------------*/

// int get_grid_decomposition(int size, int thisRank, int gridsize)
// {
//   int nbins = gridsize/size;
//   if(thisRank == size-1)
//   {
//     nbins = gridsize - nbins*(size-1);
//   }
//   
//   return nbins;
// }

/*----------------------------------------------------------------*/
/* READING GRIDS */
/*----------------------------------------------------------------*/

float *read_grid(char *filename, int gridsize, int doublePrecision, int memoryIntensive)
{
  float *thisGrid = NULL; 
  
  int zNbins = 0; //get_grid_decomposition(size, thisRank, gridsize);
  int offset = 0; //thisRank * (gridsize/size);
  
  get_grid_decomposition(gridsize, &zNbins, &offset, memoryIntensive);
  
  /* if memory efficient */
  if (memoryIntensive == 0)
  {
#ifdef MPI
    thisGrid = allocate_array_float(gridsize * gridsize * zNbins, "thisGrid");
#else
    thisGrid = allocate_array_float(gridsize * gridsize * gridsize, "thisGrid");
#endif
  }
  /* if memory inefficient but fast ? */
  else
  {
    zNbins = gridsize;
    offset = 0;
    thisGrid = allocate_array_float(gridsize * gridsize * gridsize, "thisGrid");
  }
  
  if(doublePrecision == 1)
  {
    read_grid_doubleprecision(thisGrid, gridsize, zNbins, offset, filename);
  }
  else
  {
    read_grid_singleprecision(thisGrid, gridsize, zNbins, offset, filename);
  }
    
  return thisGrid;
}

void read_grid_singleprecision(float *toThisArray, int gridsize, int zNbins, int thisOffset, char *filename)
{
#ifdef MPI
  int success;
  int resultlen;
  char msg[MPI_MAX_ERROR_STRING];

  MPI_File mpifile;
  MPI_Offset offset;
  MPI_Status status;
  
  offset = (thisOffset * gridsize * gridsize * sizeof(float));
    
  success = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &mpifile);
  if(success != MPI_SUCCESS)
  {
    MPI_Error_string(success, msg, &resultlen);
    fprintf(stderr, "MPI_File_open(): %s\n", msg);
    exit(EXIT_FAILURE);
  }
  MPI_File_read_at_all(mpifile, offset, toThisArray, zNbins * gridsize * gridsize, MPI_FLOAT, &status);
  MPI_File_close(&mpifile);
#else
  FILE *fp;
  fp = fopen(filename, "rb");
  fread(toThisArray, sizeof(float), gridsize * gridsize * gridsize, fp);
  fclose(fp);
#endif
}

void read_grid_doubleprecision(float *toThisArray, int gridsize, int zNbins, int thisOffset, char *filename)
{
  double *tmparray = allocate_array_double(zNbins*gridsize*gridsize, "tmparray");
  
#ifdef MPI
  int success;
  int resultlen;
  char msg[MPI_MAX_ERROR_STRING];

  MPI_File mpifile;
  MPI_Offset offset;
  MPI_Status status;
  
  offset = (thisOffset * gridsize * gridsize * sizeof(double));
    
  success = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,MPI_INFO_NULL, &mpifile);
  if(success != MPI_SUCCESS)
  {
    MPI_Error_string(success, msg, &resultlen);
    fprintf(stderr, "MPI_File_open(): %s\n", msg);
    exit(EXIT_FAILURE);
  }
  MPI_File_read_at_all(mpifile, offset, tmparray, zNbins * gridsize * gridsize, MPI_DOUBLE, &status);
  MPI_File_close(&mpifile);
#else
  FILE *fp;
  fp = fopen(filename, "rb");
  fread(tmparray, sizeof(double), gridsize * gridsize * gridsize, fp);
  fclose(fp);
#endif
  
  for(int i=0; i<zNbins; i++)
  {
    for(int j=0; j<gridsize; j++)
    {
      for(int k=0; k<gridsize; k++)
      {
        toThisArray[i*gridsize*gridsize+j*gridsize+k] = tmparray[i*gridsize*gridsize+j*gridsize+k];
      }
    }
  }
  free(tmparray);
}

/*----------------------------------------------------------------*/
/* COMPUTE XHI GRID FROM XHII GRID */
/*----------------------------------------------------------------*/

void get_XHI_grid(float **thisGrid, int gridsize, int memoryIntensive)
{
  int local_n0 = gridsize;
  int local_0_start = 0;
  
  get_grid_decomposition(gridsize, &local_n0, &local_0_start, memoryIntensive);
  
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<gridsize; j++)
    {
      for(int k=0; k<gridsize; k++)
      {
        (*thisGrid)[i*gridsize*gridsize+j*gridsize+k] = 1.-(*thisGrid)[i*gridsize*gridsize+j*gridsize+k];
      }
    }
  }
}

/*----------------------------------------------------------------*/
/* GETTING GRID PROPERTIES FOR GALAXIES */
/*----------------------------------------------------------------*/

int check_grid_property_present(char *propertyName)
{
  int result = 0;
  
  /* check whether a grid property is being asked for */
  if(strcmp(propertyName, "DENS") == 0)
    result = 1;
  if(strcmp(propertyName, "XHII") == 0)
    result = 2;
  if(strcmp(propertyName, "XHI") == 0)
    result = 3;
  if(strcmp(propertyName, "VELX") == 0)
    result = 4;
  if(strcmp(propertyName, "VELY") == 0)
    result = 5;
  if(strcmp(propertyName, "VELZ") == 0)
    result = 6;
  if(strcmp(propertyName, "DENSsmooth") == 0)
    result = 11;
  if(strcmp(propertyName, "XHIIsmooth") == 0)
    result = 12;
  if(strcmp(propertyName, "XHIsmooth") == 0)
    result = 13;
  if(strcmp(propertyName, "VELXsmooth") == 0)
    result = 14;
  if(strcmp(propertyName, "VELYsmooth") == 0)
    result = 15;
  if(strcmp(propertyName, "VELZsmooth") == 0)
    result = 16;
  
  return result;
}

void get_grid_values(int numEntries, double **toThisArray, double *posIndex, int gridsize, float *thisGrid, int memoryIntensive)
{
  /* figure out offset for grid and local indices */
  int local_n0 = 0, local_0_start = 0;
  get_grid_decomposition(gridsize, &local_n0, &local_0_start, memoryIntensive);

  int offset = local_0_start; //thisRank * (gridsize/size);
  int gridIndex = 0;
    
  for(int i=0; i<numEntries; i++)
  {
    gridIndex = (int)(posIndex[i]) - offset*gridsize*gridsize;
    (*toThisArray)[i] = thisGrid[gridIndex];
  }
}

double *get_grid_property(int size, int thisRank, int gridProperty, char *gridFilename, int numGal, double *posIndex, int gridsize, int doublePrecision, int memoryIntensive, int smoothedGrid, double smoothingScale)
{
  /* read in grid */
  double *thisProperty = NULL;
  
  float *thisGrid = read_grid(gridFilename, gridsize, doublePrecision, memoryIntensive);
  
  if(gridProperty == 3 || gridProperty == 13)
  {
    get_XHI_grid(&thisGrid, gridsize, memoryIntensive);
  }
  
  if(smoothedGrid == 1)
  {
    smooth_grid_property(&thisGrid, gridsize, memoryIntensive, smoothingScale);
  }
  
  /* get grid value for galaxies */
#ifdef MPI
  if(memoryIntensive == 0)
  {
    thisProperty = communicate_array_double(size, thisRank, numGal, posIndex, gridsize, posIndex, thisGrid);
  }
  else
  {
    thisProperty = allocate_array_double(numGal, "thisProperty");
    get_grid_values(numGal, &thisProperty, posIndex, gridsize, thisGrid, 1);
  }
#else
  thisProperty = allocate_array_double(numGal, "thisProperty");
  get_grid_values(numGal, &thisProperty, posIndex, gridsize, thisGrid, 1);
#endif
  
  free(thisGrid);
  
  return thisProperty;
}

double *getGridProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  char fileEnding[10];
  char *gridFilename = NULL;
  double *newProperty = NULL;
  int doublePrecision = 0;
  int smoothedGrid = 0;
  double smoothingScale = (simParam->smoothingScale / simParam->boxsize) * simParam->gridsize;
  int gridProperty = check_grid_property_present(propertyName);
  
  if(gridProperty == 1 || gridProperty == 11)
  {
    sprintf(fileEnding, "_%03d", currSnap);
    gridFilename = concat(simParam->densFilename, fileEnding);
    doublePrecision = simParam->densInputInDoublePrecision;
  }
  else if(gridProperty == 2 || gridProperty == 12)
  {
    sprintf(fileEnding, "_%02d", currSnap);
    gridFilename = concat(simParam->ionFilename, fileEnding);
    doublePrecision = simParam->ionInputInDoublePrecision;
  }
  else if(gridProperty == 3 || gridProperty == 13)
  {
    sprintf(fileEnding, "_%02d", currSnap);
    gridFilename = concat(simParam->ionFilename, fileEnding);
    doublePrecision = simParam->ionInputInDoublePrecision;
  }
  else if(gridProperty == 4 || gridProperty == 14)
  {
    sprintf(fileEnding, "_%03d", currSnap);
    gridFilename = concat(simParam->velxFilename, fileEnding);
    doublePrecision = simParam->velInputInDoublePrecision;
  }
  else if(gridProperty == 5 || gridProperty == 15)
  {
    sprintf(fileEnding, "_%03d", currSnap);
    gridFilename = concat(simParam->velyFilename, fileEnding);
    doublePrecision = simParam->velInputInDoublePrecision;
  }
  else if(gridProperty == 6 || gridProperty == 16)
  {
    sprintf(fileEnding, "_%03d", currSnap);
    gridFilename = concat(simParam->velzFilename, fileEnding);
    doublePrecision = simParam->velInputInDoublePrecision;
  }
  
  if(gridProperty > 10)
    smoothedGrid = 1;
    
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
  double *posIndex = getProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "POS", currSnap, times);
  newProperty = get_grid_property(simParam->size, simParam->thisRank, gridProperty, gridFilename, numGal, posIndex, simParam->gridsize, doublePrecision, simParam->memoryIntensive, smoothedGrid, smoothingScale);
  free(posIndex);
  free(gridFilename);
  
  return newProperty;
}
