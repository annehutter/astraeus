#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"
#include "domain.h"
#include "gal_gtree.h"

#include "cifog/confObj.h"
#include "cifog/grid.h"

#include "photion_gal.h"

/* -------------------------------------------------------- */
/* GENERAL FUNCTION TO OBTAIN THE GRID PROPERTY OF A GALAXY */
/* -------------------------------------------------------- */

double get_gridProperty_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *gridProperty)
{
  int32_t x = thisGal->pos[0] * inv_boxsize * nbins;  // assumes pos[i] to arange between 0 and 1
  int32_t y = thisGal->pos[1] * inv_boxsize * nbins;
  int32_t z = thisGal->pos[2] * inv_boxsize * nbins;
  
  if(x >= nbins)
  {
    printf("grid property: galaxy position is exactly at the box's edge: x = %e (%d)\n", thisGal->pos[0], x);
    x = nbins-1;
  }
  if(y >= nbins)
  {
    y = nbins-1;
    printf("grid property: galaxy position is exactly at the box's edge: y = %e (%d)\n", thisGal->pos[1], y);
  }
  if(z >= nbins)
  {
    z = nbins-1;
    printf("grid property: galaxy position is exactly at the box's edge: z = %e (%d)\n", thisGal->pos[2], z);
  }
 
  int32_t index = z*nbins*nbins + y*nbins + x;
  
  if(gridProperty != NULL)
    return creal(gridProperty[index]);
  else
    return 0.;
}

/* --------------------------------- */
/* GENERAL FUNCTIONS TO UPDATE FIELD */
/* --------------------------------- */

void update_field(fftw_complex **field, grid_t *thisGrid, double redshift, double XHII_threshold, char *fieldName, int memoryIntensive)
{ 
  update_grid_field(thisGrid, redshift, XHII_threshold, fieldName);
  
  if(memoryIntensive == 1)
  {
    update_local_field(field, thisGrid, XHII_threshold, fieldName);
  }
}

void update_grid_field(grid_t *thisGrid, double redshift, double XHII_threshold, char *fieldName)
{
  int32_t local_n0 = thisGrid->local_n0;
  int32_t nbins = thisGrid->nbins;
  fftw_complex *thisGridField = NULL;
  int32_t fieldIdentifier = -1;
  
  if(strcmp(fieldName, "ZREION") == 0)
  {
    thisGridField = thisGrid->zreion;
    fieldIdentifier = 0;
  }
  else if (strcmp(fieldName, "PHOTHI") == 0)
  {
    thisGridField = thisGrid->photHI_zreion;
    fieldIdentifier = 1;
  }
  else
  {
    fprintf(stderr, "Not valid option for field!\n");
    exit(EXIT_FAILURE);
  }
  
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        if(creal(thisGrid->zreion[i*nbins*nbins+j*nbins+k]) <= 0. && creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k]) > XHII_threshold)
        {
          if(fieldIdentifier == 0)
            thisGridField[i*nbins*nbins+j*nbins+k] = redshift + 0.*I;
          if(fieldIdentifier == 1)
            thisGridField[i*nbins*nbins+j*nbins+k] = thisGrid->photHI[i*nbins*nbins+j*nbins+k];
        }
      }
    }
  }
}

void update_local_field(fftw_complex **field, grid_t *thisGrid, double XHII_threshold, char *fieldName)
{
  int32_t local_n0 = thisGrid->local_n0;
  int32_t nbins = thisGrid->nbins;
  fftw_complex *thisGridField = NULL;
  
  if(strcmp(fieldName, "ZREION") == 0)
  {
    thisGridField = thisGrid->zreion;
  }
  else if (strcmp(fieldName, "PHOTHI") == 0)
  {
    thisGridField = thisGrid->photHI_zreion;
  }
  else
  {
    fprintf(stderr, "Not valid option for field!\n");
    exit(EXIT_FAILURE);
  }
  
  if(*field == NULL)
  {
    *field = fftw_malloc(sizeof(fftw_complex) * nbins*nbins*nbins);
    initialize_grid(*field, thisGrid->nbins, thisGrid->local_n0, -1.);
  }

#ifdef MPI    
  MPI_Allgather(thisGridField, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, *field, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
#else
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        if(creal(thisGrid->zreion[i*nbins*nbins+j*nbins+k]) <= 0. && creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k]) > XHII_threshold)
        {
          (*field)[i*nbins*nbins+j*nbins+k] = thisGridField[i*nbins*nbins+j*nbins+k];
        }
      }
    }
  }
#endif
}

void write_grid_to_file_double_test(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename, int thisRank, int memory_intensive)
{
    double *tmparray;
#ifdef MPI
    MPI_File mpifile;
    MPI_Offset offset;
    MPI_Status status;
#endif
    
    tmparray = (double*)malloc(sizeof(double)*local_n0*nbins*nbins);
    if(tmparray == NULL)
    {
        fprintf(stderr, "tmparray in write_grid_to_file_double (grid.c) could not be allocated\n");
        exit(EXIT_FAILURE);
    }
    
    for(int i=0; i<local_n0; i++)
    {
        for(int j=0; j<nbins; j++)
        {
            for(int k=0; k<nbins; k++)
            {
                tmparray[i*nbins*nbins+j*nbins+k] = (double)creal(thisArray[i*nbins*nbins+j*nbins+k]);
            }
        }
    }
    
#ifdef MPI
    if(memory_intensive == 0)
    {
        offset = (local_0_start*nbins*nbins*sizeof(double));
        
        MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifile);
        if (mpifile == MPI_FILE_NULL)
        {
            fprintf(stderr, "Could not open file %s\n", filename);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
        }
        MPI_File_write_at_all(mpifile, offset, tmparray, local_n0*nbins*nbins, MPI_DOUBLE, &status);
        MPI_File_close(&mpifile);
    }
#endif
    
    if(thisRank == 0 && memory_intensive == 1)
    {
        FILE * fp;
        
        fp = fopen(filename, "wb");
        if (fp == NULL)
        {
            fprintf(stderr, "Could not open file %s\n", filename);
            exit(EXIT_FAILURE);
        } 
        fwrite(tmparray, sizeof(double), nbins*nbins*nbins, fp);
        fclose(fp);
    }
    
    free(tmparray);
}