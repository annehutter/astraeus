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
#include "communicate_mpi.h"
#include "get_grid_properties.h"

/*----------------------------------------------------------------*/
/* DO DOMAIN DECOMPOSITION FOR GRIDS */
/*----------------------------------------------------------------*/

int get_grid_decomposition(int size, int thisRank, int gridsize)
{
  int nbins = gridsize/size;
  if(thisRank == size-1)
  {
    nbins = gridsize - nbins*(size-1);
  }
  
  return nbins;
}

/*----------------------------------------------------------------*/
/* READING GRIDS */
/*----------------------------------------------------------------*/

float *read_grid(int size, int thisRank, char *filename, int gridsize, int doublePrecision, int memoryIntensive)
{
  float *thisGrid = NULL; 
  
  int zNbins = get_grid_decomposition(size, thisRank, gridsize);
  int offset = thisRank * (gridsize/size);
  
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
  
  offset = (thisOffset * gridsize * gridsize * sizeof(double));
    
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
/* GETTING GRID PROPERTIES FOR GALAXIES */
/*----------------------------------------------------------------*/

int check_grid_property_present(char *propertyName)
{
  int result = 0;
  
  /* check whether a grid property is being asked for */
  if(strcmp(propertyName, "DENS") == 0)
    result = 1;
  if(strcmp(propertyName, "ION") == 0)
    result = 2;
  
  return result;
}

float *get_grid_property(int size, int thisRank, char *gridFilename, int numGal, float *posIndex, int gridsize, int doublePrecision, int memoryIntensive)
{
  /* read in grid */
//   int doublePrecision = 0;
//   int memoryIntensive = 0;
  float *thisProperty = NULL;
  
  float *thisGrid = read_grid(size, thisRank, gridFilename, gridsize, doublePrecision, memoryIntensive);
  
  /* get grid value for galaxies */
  if(memoryIntensive == 0)
  {
    thisProperty = communicate_array_float(size, thisRank, numGal, posIndex, gridsize, posIndex, thisGrid);
  }
  else
  {
    thisProperty = allocate_array_float(numGal, "thisProperty");
    get_grid_values(1, 0, numGal, &thisProperty, posIndex, gridsize, thisGrid);
  }
  
  free(thisGrid);
  
  return thisProperty;
}

void get_grid_values(int size, int thisRank, int numEntries, float **toThisArray, float *posIndex, int gridsize, float *thisGrid)
{
  /* figure out offset for grid and local indices */
//   int zNbins = get_grid_decomposition(size, thisRank, gridsize);
  int offset = thisRank * (gridsize/size);
  int gridIndex = 0;
    
  for(int i=0; i<numEntries; i++)
  {
    gridIndex = (int)(posIndex[i]) - offset*gridsize*gridsize;
    (*toThisArray)[i] = thisGrid[gridIndex];
  }
}
