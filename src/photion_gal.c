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

void update_photHI_field(fftw_complex **photHI, grid_t *thisGrid)
{
  int32_t local_n0 = thisGrid->local_n0;
  int32_t nbins = thisGrid->nbins;
  
  if(*photHI == NULL)
  {
    *photHI = fftw_malloc(sizeof(fftw_complex) * nbins*nbins*nbins);
  }

#ifdef MPI  
  MPI_Allgather(thisGrid->photHI, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, *photHI, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
#else
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        (*photHI)[i*nbins*nbins+j*nbins+k] = thisGrid->photHI[i*nbins*nbins+j*nbins+k];
      }
    }
  }
#endif
}

void update_photHI_zreion_field(fftw_complex **photHI_zreion, grid_t *thisGrid, double redshift, double XHII_threshold)
{
  int32_t local_0_start = thisGrid->local_0_start;
  int32_t local_n0 = thisGrid->local_n0;
  int32_t nbins = thisGrid->nbins;
  
  if(*photHI_zreion == NULL)
  {
    *photHI_zreion = fftw_malloc(sizeof(fftw_complex) * nbins*nbins*nbins);
    initialize_grid(*photHI_zreion, thisGrid->nbins, thisGrid->local_n0, -1.);
  }

#ifdef MPI 
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        if(creal(thisGrid->zreion[i*nbins*nbins+j*nbins+k]) <= 0. && creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k]) > XHII_threshold)
          thisGrid->photHI_zreion[i*nbins*nbins+j*nbins+k] = thisGrid->photHI[i*nbins*nbins+j*nbins+k];
      }
    }
  }
  
  MPI_Allgather(thisGrid->photHI_zreion, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, *photHI_zreion, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
#else
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        if(creal(thisGrid->zreion[i*nbins*nbins+j*nbins+k]) <= 0. && creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k]) > XHII_threshold)
          (*photHI_zreion)[i*nbins*nbins+j*nbins+k] = thisGrid->photHI_zreion[i*nbins*nbins+j*nbins+k];
      }
    }
  }
#endif
}

double get_photHI_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *photHI)
{ 
  int32_t x = thisGal->pos[0] * inv_boxsize * nbins;  // assumes pos[i] to arange between 0 and 1
  int32_t y = thisGal->pos[1] * inv_boxsize * nbins;
  int32_t z = thisGal->pos[2] * inv_boxsize * nbins;
  
  if(x >= nbins)
  {
    printf("photHI: galaxy position is exactly at the box's edge: x = %e (%d)\n", thisGal->pos[0], x);
    x = nbins-1;
  }
  if(y >= nbins)
  {
    y = nbins-1;
    printf("photHI: galaxy position is exactly at the box's edge: y = %e (%d)\n", thisGal->pos[1], y);
  }
  if(z >= nbins)
  {
    z = nbins-1;
    printf("photHI: galaxy position is exactly at the box's edge: z = %e (%d)\n", thisGal->pos[2], z);
  }
 
  int32_t index = z*nbins*nbins + y*nbins + x;
  
  if(photHI != NULL)
    return creal(photHI[index]);
  else
    return 0.;
}

/* ---------------------------------------- */


void update_XHII_field(fftw_complex **XHII, grid_t *thisGrid)
{
  int32_t local_n0 = thisGrid->local_n0;
  int32_t nbins = thisGrid->nbins;
  
  if(*XHII == NULL)
  {
    *XHII = fftw_malloc(sizeof(fftw_complex) * nbins*nbins*nbins);
  }

#ifdef MPI  
  MPI_Allgather(thisGrid->XHII, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, *XHII, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
#else
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        (*XHII)[i*nbins*nbins+j*nbins+k] = thisGrid->XHII[i*nbins*nbins+j*nbins+k];
      }
    }
  }
#endif
}

double get_XHII_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *XHII)
{ 
  int32_t x = thisGal->pos[0] * inv_boxsize * nbins;  // assumes pos[i] to arange between 0 and 1
  int32_t y = thisGal->pos[1] * inv_boxsize * nbins;
  int32_t z = thisGal->pos[2] * inv_boxsize * nbins;

  if(x >= nbins)
  {
    printf("XHII: galaxy position is exactly at the box's edge: x = %e (%d)\n", thisGal->pos[0], x);
    x = nbins-1;
  }
  if(y >= nbins)
  {
    y = nbins-1;
    printf("XHII: galaxy position is exactly at the box's edge: y = %e (%d)\n", thisGal->pos[1], y);
  }
  if(z >= nbins)
  {
    z = nbins-1;
    printf("XHII: galaxy position is exactly at the box's edge: z = %e (%d)\n", thisGal->pos[2], z);
  }
  
  int32_t index = z*nbins*nbins + y*nbins + x;
  
  if(XHII != NULL)
    return creal(XHII[index]);
  else
    return 0.;
}

/* ---------------------------------------- */

void update_zreion_field(fftw_complex **zreion, grid_t *thisGrid, double redshift, double XHII_threshold)
{
  int32_t local_0_start = thisGrid->local_0_start;
  int32_t local_n0 = thisGrid->local_n0;
  int32_t nbins = thisGrid->nbins;
  
  if(*zreion == NULL)
  {
    *zreion = fftw_malloc(sizeof(fftw_complex) * nbins*nbins*nbins);
    initialize_grid(*zreion, thisGrid->nbins, thisGrid->local_n0, -1.);
  }

#ifdef MPI 
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        if(creal(thisGrid->zreion[i*nbins*nbins+j*nbins+k]) <= 0. && creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k]) > XHII_threshold)
          thisGrid->zreion[i*nbins*nbins+j*nbins+k] = redshift + 0.*I;
      }
    }
  }
  
  MPI_Allgather(thisGrid->zreion, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, *zreion, local_n0*nbins*nbins, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
#else
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        if(creal(thisGrid->zreion[i*nbins*nbins+j*nbins+k]) <= 0. && creal(thisGrid->XHII[i*nbins*nbins+j*nbins+k]) > XHII_threshold)
          (*zreion)[i*nbins*nbins+j*nbins+k] = redshift + 0*I;
      }
    }
  }
#endif
}

double get_zreionField_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *zreion)
{ 
  int32_t x = thisGal->pos[0] * inv_boxsize * nbins;  // assumes pos[i] to arange between 0 and 1
  int32_t y = thisGal->pos[1] * inv_boxsize * nbins;
  int32_t z = thisGal->pos[2] * inv_boxsize * nbins;

  if(x >= nbins)
  {
    printf("zreion: galaxy position is exactly at the box's edge: x = %e (%d)\n", thisGal->pos[0], x);
    x = nbins-1;
  }
  if(y >= nbins)
  {
    y = nbins-1;
    printf("zreion: galaxy position is exactly at the box's edge: y = %e (%d)\n", thisGal->pos[1], y);
  }
  if(z >= nbins)
  {
    z = nbins-1;
    printf("zreion: galaxy position is exactly at the box's edge: z = %e (%d)\n", thisGal->pos[2], z);
  }
  
  int32_t index = z*nbins*nbins + y*nbins + x;
  
  if(zreion != NULL)
    return creal(zreion[index]);
  else
    return 0.;
}

double get_zreion_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *XHII, double XHII_threshold)
{
  double XHII_gal = get_XHII_gal(thisGal, nbins, inv_boxsize, XHII);
  double zreion = 0;
  
  if(thisGal->zreion == 0 && XHII_gal > XHII_threshold)
  {
    zreion = 1./thisGal->scalefactor - 1.;
  }

  return zreion;
}


void write_grid_to_file_double_test(fftw_complex *thisArray, int nbins, int local_n0, char *filename)
{
    double *tmparray;
    
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
    
    FILE * fp;
    
    fp = fopen(filename, "wb");
    if (fp == NULL)
    {
        fprintf(stderr, "Could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    } 
    fwrite(tmparray, sizeof(double), nbins*nbins*nbins, fp);
    fclose(fp);
    
    free(tmparray);
}