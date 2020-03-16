#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"
#include "nion.h"

nion_t *initNion(int32_t numGal)
{
  nion_t *thisNion = malloc(sizeof(nion_t));
  
  if(thisNion == NULL)
  {
    fprintf(stderr, "Could not allocate nion_t\n");
    exit(EXIT_FAILURE);
  }
  
  thisNion->numGal = numGal;
  thisNion->numGalWritten = 0;
  thisNion->nion = allocate_array_double(numGal, "nion");
  thisNion->pos = allocate_array_int32_t(numGal, "pos");
  thisNion->rank = allocate_array_int32_t(numGal, "rank");

  return thisNion;
}

void reallocNion(nion_t **thisNion, int32_t numGal)
{  
  (*thisNion)->numGal = numGal;
  (*thisNion)->nion = realloc((*thisNion)->nion, sizeof(double) * numGal);
  (*thisNion)->pos = realloc((*thisNion)->pos, sizeof(int32_t) * numGal);
  (*thisNion)->rank = realloc((*thisNion)->rank, sizeof(int32_t) * numGal);
}

void deallocate_nion(nion_t *thisNion)
{
  if(thisNion->nion != NULL) free(thisNion->nion);
  if(thisNion->pos != NULL) free(thisNion->pos);
  if(thisNion->rank != NULL) free(thisNion->rank);

  free(thisNion);
}

double get_mean_nion(nion_t *thisNion)
{
  double sum = 0.;
  
  for(int i=0; i<thisNion->numGalWritten; i++)
  {
    sum += 1.e-50*thisNion->nion[i];
  }
  
#ifdef MPI
  double sum_all = 0.;
  MPI_Allreduce(&sum, &sum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  sum = sum_all;
#endif
  
  return sum;
}