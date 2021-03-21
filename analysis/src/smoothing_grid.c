#include <stdio.h>
#include <stdlib.h>
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

#include "utils.h"
#include "dconfObj.h"
#include "smoothing_grid.h"

#define SQR(X) ((X) * (X))

void smooth_grid_property(float **thisGrid, int nbins, int memoryIntensive, float smooth_scale)
{
  ptrdiff_t local_n0 = nbins;
  ptrdiff_t local_0_start = 0;
  fftw_complex *thisInput = NULL;
  fftw_complex *thisSmoothedGrid = NULL;
    
#ifdef MPI
  ptrdiff_t alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
  
  if(memoryIntensive == 0)
  {
    thisInput = fftw_alloc_complex(alloc_local);
    thisSmoothedGrid = fftw_alloc_complex(alloc_local);
  }
  else
  {
    local_n0 = nbins;
    local_0_start = 0;
    thisInput = fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    thisSmoothedGrid = fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
  }
#else  
  thisInput = fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
  thisSmoothedGrid = fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
#endif
  
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        thisInput[i*nbins*nbins+j*nbins+k] = (*thisGrid)[i*nbins*nbins+j*nbins+k] + 0.*I;
        thisSmoothedGrid[i*nbins*nbins+j*nbins+k] = 0. + 0.*I;
      }
    }
  }
  
  fftw_complex *filter = construct_tophat_filter(nbins, local_0_start, local_n0, smooth_scale, memoryIntensive);
  
  convolve_fft(nbins, local_0_start, local_n0, filter, &thisSmoothedGrid, &thisInput, memoryIntensive);

  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        (*thisGrid)[i*nbins*nbins+j*nbins+k] = creal(thisSmoothedGrid[i*nbins*nbins+j*nbins+k]);
      }
    }
  }
  
  fftw_free(filter);
  fftw_free(thisInput);
  fftw_free(thisSmoothedGrid);
}

void get_grid_decomposition(int nbins, int *local_n0, int *local_0_start, int memoryIntensive)
{
  *local_n0 = nbins;
  *local_0_start = 0;
  
#ifdef MPI
  ptrdiff_t local_n0_local, local_0_start_local;
  
  if(memoryIntensive == 0)
  {
    fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0_local, &local_0_start_local);
    
    *local_n0 = local_n0_local;
    *local_0_start = local_0_start_local;
  }
#endif
}

fftw_complex *construct_tophat_filter(int nbins, ptrdiff_t local_0_start, ptrdiff_t local_n0, float smooth_scale, int memoryIntensive)
/* nbins: gridsize of the fftw field */
/* local_0_start: offset of fftw field along z axis */
/* local_n0: size of fftw field along z axis */
/* smooth_scale: smoothing scale in units of bins */
{
  double normCoeff = 0.;
  const double half_nbins = nbins*0.5;
  const double sq_smooth_scale = 0.25*SQR(smooth_scale);
  fftw_complex *filter = NULL;
  
#ifdef MPI
  ptrdiff_t alloc_local, local_n0_local, local_0_start_local;

  if(memoryIntensive == 0)
  {
    alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0_local, &local_0_start_local);
    assert(local_n0 == local_n0_local);
    assert(local_0_start == local_0_start_local);
    
    filter = fftw_alloc_complex(alloc_local);
  }
  else
  {
    filter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
  }
#else
  filter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
#endif

  for(int i=0; i<nbins; i++)
  {      
    const double i_expr = half_nbins-fabs(i - half_nbins);        
    const double sq_i_expr = SQR(i_expr);
    for(int j=0; j<nbins; j++)
    {
      const double j_expr = half_nbins - fabs(j - half_nbins);
      const double sq_j_expr = SQR(j_expr);
      for(int k=0; k<nbins; k++)
      {
        const double k_expr = half_nbins - fabs(k - half_nbins);
        const double sq_k_expr = SQR(k_expr);
        
        const double expr = sq_i_expr + sq_j_expr + sq_k_expr - sq_smooth_scale;
        normCoeff += (expr<=0.) ? 1.:0.;
      }
    }
  }
  normCoeff = 1./normCoeff;
  
  for(int i=0; i<local_n0; i++)
  {
    const double i_expr = half_nbins - fabs(i + local_0_start - half_nbins);        
    const double sq_i_expr = SQR(i_expr);
    for(int j=0; j<nbins; j++)
    {
      const double j_expr = half_nbins - fabs(j - half_nbins);
      const double sq_j_expr = SQR(j_expr);
      for(int k=0; k<nbins; k++)
      {            
        const double k_expr = half_nbins - fabs(k - half_nbins);
        const double sq_k_expr = SQR(k_expr);
      
        const double expr = sq_i_expr + sq_j_expr + sq_k_expr - sq_smooth_scale;
        
        filter[i*nbins*nbins+j*nbins+k] = (expr<=0.) ? normCoeff+0.*I : 0.+0.*I;
      }
    }
  }
  
  return filter;
}

void convolve_fft(int nbins, ptrdiff_t local_0_start, ptrdiff_t local_n0, fftw_complex *filter, fftw_complex **output, fftw_complex **input, int memoryIntensive)
{
  double factor;
#ifdef MPI
  ptrdiff_t alloc_local, local_n0_local, local_0_start_local;
#else
  ptrdiff_t local_n0_local;
#endif
  
  fftw_complex *input_ft, *filter_ft;
  fftw_plan plan_input, plan_filter, plan_back;
  
  local_n0_local = local_n0;
  
#ifdef MPI
  if(memoryIntensive == 0)
  {
    local_0_start_local = local_0_start;
    
    alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0_local, &local_0_start_local);
    assert(local_0_start == local_0_start);
    assert(local_n0 == local_n0);
    
    input_ft = fftw_alloc_complex(alloc_local);
    filter_ft = fftw_alloc_complex(alloc_local);
    
    plan_input = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, *input, input_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE); //FFTW_MPI_TRANSPOSED_OUT); //FFTW_MPI_TRANSPOSED_OUT
    plan_filter = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE); //FFTW_MPI_TRANSPOSED_OUT);
  }
  else
  {
    input_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    if(input_ft == NULL)
    {
      fprintf(stderr, "input_ft in convolve_fft (filtering.c) could not be allocated\n");
      exit(EXIT_FAILURE);
    }
    filter_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
    if(filter_ft == NULL)
    {
      fprintf(stderr, "filter_ft in convolve_fft (filtering.c) could not be allocated\n");
      exit(EXIT_FAILURE);
    }
    
    plan_input = fftw_plan_dft_3d(nbins, nbins, nbins, *input, input_ft, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_filter = fftw_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, FFTW_FORWARD, FFTW_ESTIMATE);
  }
#else 
  input_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
  if(input_ft == NULL)
  {
    fprintf(stderr, "input_ft in convolve_fft (filtering.c) could not be allocated\n");
    exit(EXIT_FAILURE);
  }
  filter_ft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
  if(filter_ft == NULL)
  {
    fprintf(stderr, "filter_ft in convolve_fft (filtering.c) could not be allocated\n");
    exit(EXIT_FAILURE);
  }
  
  plan_input = fftw_plan_dft_3d(nbins, nbins, nbins, *input, input_ft, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_filter = fftw_plan_dft_3d(nbins, nbins, nbins, filter, filter_ft, FFTW_FORWARD, FFTW_ESTIMATE);
#endif
  
  fftw_execute(plan_input);
  fftw_execute(plan_filter);
  
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        input_ft[i*nbins*nbins+j*nbins+k] = input_ft[i*nbins*nbins+j*nbins+k]*filter_ft[i*nbins*nbins+j*nbins+k];
      }
    }
  }
  
#ifdef MPI
  if(memoryIntensive == 0)
  {
    plan_back = fftw_mpi_plan_dft_3d(nbins, nbins, nbins, input_ft, *output, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); //FFTW_MPI_TRANSPOSED_IN);
  }
  else
  {
    plan_back = fftw_plan_dft_3d(nbins, nbins, nbins, input_ft, *output, FFTW_BACKWARD, FFTW_ESTIMATE);
  }
#else
  plan_back = fftw_plan_dft_3d(nbins, nbins, nbins, input_ft, *output, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
  fftw_execute(plan_back);
  
  factor = 1./(nbins*nbins*nbins);
  
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        (*output)[i*nbins*nbins+j*nbins+k] = factor*(*output)[i*nbins*nbins+j*nbins+k];
      }
    }
  }       
  
  fftw_destroy_plan(plan_input);
  fftw_destroy_plan(plan_filter);
  fftw_destroy_plan(plan_back);
  
  fftw_free(input_ft);
  fftw_free(filter_ft);
}