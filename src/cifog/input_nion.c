#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <assert.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "confObj.h"
#include "grid.h"
#include "input_nion.h"

/* read in / update sources or nion -------------------------------------------------------------*/
void read_update_nion(confObj_t simParam, fftw_complex *nion, grid_t *thisGrid)
{
//   if(thisGrid->nion != NULL) free(thisGrid->nion);
//   thisGrid->nion = nion;
  
  int nbins = thisGrid->nbins;
  int local_n0 = thisGrid->local_n0;
  
  double sum = 0.;
  double sum_all = 0.;
  
  for(int i=0; i<local_n0; i++)
  {
      for(int j=0; j<nbins; j++)
      {
          for(int k=0; k<nbins; k++)
          {
              thisGrid->nion[i*nbins*nbins+j*nbins+k] = 1.e50*creal(nion[i*nbins*nbins+j*nbins+k]) + 0.*I;
              sum += 1.e-50*creal(thisGrid->nion[i*nbins*nbins+j*nbins+k]);
          }
      }
  }
  
#ifdef MPI
  MPI_Allreduce(&sum, &sum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  sum_all = sum;
#endif
  if(thisGrid->local_0_start == 0) printf("\n********\nCIFOG: sum(nion) = %e\n********\n", sum_all);
}

void read_update_nion_HeI(confObj_t simParam, fftw_complex *nion_HeI, grid_t *thisGrid)
{
  if(thisGrid->nion_HeI != NULL) free(thisGrid->nion_HeI);
  thisGrid->nion_HeI = nion_HeI;
}

void read_update_nion_HeII(confObj_t simParam, fftw_complex *nion_HeII, grid_t *thisGrid)
{
  if(thisGrid->nion_HeII != NULL) free(thisGrid->nion_HeII);
  thisGrid->nion_HeII = nion_HeII;
}