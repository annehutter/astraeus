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
#include "dconfObj.h"
#include "gal_gtree.h"
#include "nion.h"
#include "domain.h"
#include "ion_emissivity.h"

#include "sort.h"
#include "nion_gal.h"

/*-----------------------------------------------------*/
/* MAPPING GALAXIES ONTO GRID                          */
/*-----------------------------------------------------*/

/* GENERATE FFTW ARRAY FOR NION */
fftw_complex *map_galnion_to_grid(nion_t *thisNion, domain_t *thisDomain, int memoryIntensive)
{
#ifdef MPI
  ptrdiff_t alloc_local, local_n0, local_0_start;
#else
  ptrdiff_t local_n0, local_0_start;
#endif
  
  fftw_complex *nion = NULL;
  int nbins = thisDomain->nbins;
  local_n0 = thisDomain->local_n0;
  local_0_start = thisDomain->local_0_start;
  
#ifdef MPI  
  fftw_mpi_init();
  alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);

  assert((int)local_n0 == thisDomain->local_n0);
  assert((int)local_0_start == thisDomain->local_0_start);
  
  nion = fftw_alloc_complex(alloc_local);
#else
  nion = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nbins*nbins*nbins);
#endif
    
  if(memoryIntensive == 1)
  {
    map_nion_to_pos_on_grid_memintensive(thisNion, nion, nbins, local_n0, local_0_start, thisDomain->thisRank);
  }
  else
  {
#ifdef MPI
    distribute_nion_to_processors(thisNion, thisDomain);
#endif
    map_nion_to_pos_on_grid(thisNion, nion, nbins, local_n0, local_0_start, thisDomain->thisRank);
  }

#ifdef MPI
  fftw_mpi_cleanup();
#endif

  return nion;
}

/* MAPPING FUNCTION TO GRID */
void map_nion_to_pos_on_grid(nion_t *thisNion, fftw_complex *nion, int nbins, ptrdiff_t local_n0, ptrdiff_t local_0_start, int thisRank)
{
  int32_t numGal = thisNion->numGal;
  double *nionGal = thisNion->nion;
  int32_t *posGal = thisNion->pos;
  
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        nion[i*nbins*nbins+j*nbins+k] = 0. + 0.*I;
      }
    }
  }
  
  for(int gal=0; gal<numGal; gal++)
  {
    if(posGal[gal] - local_0_start*nbins*nbins < 0 || posGal[gal] - local_0_start*nbins*nbins > local_n0*nbins*nbins)
    {
      printf("rank %d: gal = %d: %d < pos = %d < %d\n", thisRank, gal, local_0_start*nbins*nbins, posGal[gal], (local_0_start + local_n0)*nbins*nbins);
      fprintf(stderr, "Galaxy position is not within domain. Check your mapping!\n");
      exit(EXIT_FAILURE);
    }

    nion[posGal[gal] - local_0_start*nbins*nbins] += nionGal[gal] + 0.*I;
  }
}

/* MEMORY INTENSIVE MAPPING FUNCTION TO GRID */
void map_nion_to_pos_on_grid_memintensive(nion_t *thisNion, fftw_complex *nion, int nbins, ptrdiff_t local_n0, ptrdiff_t local_0_start, int thisRank)
{   
  int32_t numGal = thisNion->numGal;
  double *nionGal = thisNion->nion;
  int32_t *posGal = thisNion->pos;
  
  double *nionGlobal = allocate_array_double(nbins*nbins*nbins, "nion grid");
#ifdef MPI
  double *nionLocal = allocate_array_double(nbins*nbins*nbins, "nion grid");
#endif

#ifdef MPI
  for(int gal=0; gal<numGal; gal++)
    nionLocal[posGal[gal]] += nionGal[gal];

  MPI_Allreduce(nionLocal, nionGlobal, nbins*nbins*nbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  free(nionLocal);
#else
  for(int gal=0; gal<numGal; gal++)
    nionGlobal[posGal[gal]] += nionGal[gal];
#endif
  
  for(int i=0; i<local_n0; i++)
  {
    for(int j=0; j<nbins; j++)
    {
      for(int k=0; k<nbins; k++)
      {
        nion[i*nbins*nbins+j*nbins+k] = nionGlobal[(i+local_0_start)*nbins*nbins+j*nbins+k] + 0.*I;
      }
    }
  }
  
  if(nionGlobal != NULL) free(nionGlobal);
}

/*=======================================================================================*/

/*-----------------------------------------------------*/
/* COMMUNICATION FOR NION OF GALAXIES TO BE MAPPED     */
/*-----------------------------------------------------*/

/* MAIN ROUTINE: SEND & RECIVE GALAXIES TO PROCESSORS ACCORDING TO THEIR LOCATION */
#ifdef MPI
void distribute_nion_to_processors(nion_t *thisNion, domain_t *thisDomain)
{
  int32_t size = thisDomain->size;
  
  double *toSendNion = NULL;
  int32_t *toSendPos = NULL;
  double *toRecvNion = NULL;
  int32_t *toRecvPos = NULL;
    
  /* SORTING GALAXIES ACCORDING TO RANKS WHERE THEY ARE SENT TO */
  sort_nion(thisNion, thisDomain->thisRank);
 
  /* NUMBER OF GALAXIES ON RANK TO SEND TO OTHER RANKS */
  int32_t *numGalToRanks = get_nion_numbers(thisNion, size);
   
  /* NUMBER OF GALAXIES THAT ARE RECEIVED BY THE ALL RANKS ON THIS RANK */
  int32_t *numGalOnThisRank = allocate_array_int32_t(size, "numGalOnThisRank");
  
//   for(int i=0; i<size; i++)
//     printf("rank %d: sending %d galaxies to rank %d\n", thisDomain->thisRank, numGalToRanks[i], i);
  
  /* copy NION & POS to temporary arrays */
  toSendNion = copy_nion(thisNion);
  toSendPos = copy_pos(thisNion);
   
  /* SENDING & RECEIVING NION & POS TO RELEVANT PROCESSORS */
  send_recv_nion_pos(thisDomain, numGalToRanks, &numGalOnThisRank, toSendNion, toSendPos, &toRecvNion, &toRecvPos);
  
  free(toSendNion);
  free(toSendPos);
  
  /* REPLACE NION & POS BY GALAXIES ON DOMAIN */
  if(thisNion->nion != NULL)
  {
    free(thisNion->nion);
    thisNion->nion = toRecvNion;
  }
  if(thisNion->pos != NULL)
  {
    free(thisNion->pos);
    thisNion->pos = toRecvPos;
  }

  int32_t sum = 0;
  for(int i=0; i<size; i++)
  {
    sum += numGalOnThisRank[i];
//     printf("rank %d: received %d galaxies from rank %d\n", thisDomain->thisRank, numGalOnThisRank[i], i);
  }
  thisNion->numGal = sum;
  thisNion->numGalWritten = sum;
  
  double sumDouble = get_mean_nion(thisNion);
  if(thisDomain->thisRank == 0) printf("nion = %e\n", sumDouble);
  
//   printf("rank %d: total number of galaxies = %d and %d\n", thisDomain->thisRank, thisNion->numGal, thisNion->numGalWritten);
  
  free(numGalToRanks);
  free(numGalOnThisRank);
}
#endif

/* SORT ARRAY ACCORDING TO RANK WHERE TO SEND IT */
void sort_nion(nion_t *thisNion, int thisRank)
{
  int32_t numGal = thisNion->numGal;
  int32_t *indexArray = allocate_array_int32_t(numGal, "indexArray");
  int32_t *posArray = allocate_array_int32_t(numGal, "posArray");
  double *nionArray = allocate_array_double(numGal, "nionArray");
  
  for(int gal=0; gal<numGal; gal++) indexArray[gal] = gal;

//   if(thisRank == 0)
//   {
//     printf("before: ");
//     for(int gal=0; gal<numGal; gal++) printf("%d\t", indexArray[gal]);
//     printf("\n");
//   }
  
  quicksort(thisNion->rank, 0, numGal-1, indexArray);

//   if(thisRank == 0)
//   {
//     printf("after: ");
//     for(int gal=0; gal<numGal; gal++) printf("%d\t", indexArray[gal]);
//     printf("\n");
//   }
  
  for(int gal=0; gal<numGal; gal++)
  {
    posArray[gal] = thisNion->pos[indexArray[gal]];
    nionArray[gal] = thisNion->nion[indexArray[gal]];
  }
  
  free(indexArray);
  free(thisNion->pos);
  free(thisNion->nion);
  
  thisNion->pos = posArray;
  thisNion->nion = nionArray;
}

/* GET HOW MANY GALAXIES HAVE TO BE SENT TO EACH RANK (works only on a sorted array!) */
int32_t *get_nion_numbers(nion_t *thisNion, int32_t size)
{
  int32_t numGal = thisNion->numGal;
  int32_t *numGalToRanks = allocate_array_int32_t(size, "numGalToRanks");
  
  for(int gal=0; gal<numGal; gal++)
  {
    numGalToRanks[thisNion->rank[gal]]++;
  }
    
  return numGalToRanks;
}

/* GENERATE A COPY OF NION IN THISNION */
double *copy_nion(nion_t *thisNion)
{
  int32_t numGal = thisNion->numGal;
  double *nion = allocate_array_double(numGal, "nion");
  
  for(int gal=0; gal<numGal; gal++)
  {
    nion[gal] = thisNion->nion[gal];
  }
    
  return nion;
}

/* GENERATE A COPY OF POS IN THISNION */
int32_t *copy_pos(nion_t *thisNion)
{
  int32_t numGal = thisNion->numGal;
  int32_t *pos = allocate_array_int32_t(numGal, "pos");
  
  for(int gal=0; gal<numGal; gal++)
  {
    pos[gal] = thisNion->pos[gal];
  }
  
  return pos;
}

/* SEND & RECEIVE NION & POS TO PROCESSORS ACCORDING TO THEIR LOCATION */
#ifdef MPI
void send_recv_nion_pos(domain_t *thisDomain, int32_t *numGalToRanks, int32_t **numGalOnThisRank, double *toSendNion, int32_t *toSendPos, double **toRecvNion, int32_t **toRecvPos)
{
  MPI_Status status;
  int size = thisDomain->size;
  int thisRank = thisDomain->thisRank;
  
  int32_t head_offset = 0;
  
  int32_t *recvNumGalOnRank = *numGalOnThisRank;
  int32_t counterRecvNumGal = 0;
  
  /* head is the sending rank */
  for(int head=0; head<size; head++)
  {    
    head_offset = 0;
    for(int destRank=0; destRank<size; destRank++)
    {
      if((thisRank == head) && (thisRank != destRank))
      {
        /* send stuff from head to destRank */
        MPI_Send(&numGalToRanks[destRank], 1, MPI_INT, destRank, 100, MPI_COMM_WORLD);
        
        if(numGalToRanks[destRank] > 0)
        {            
          MPI_Send(&toSendNion[head_offset], numGalToRanks[destRank], MPI_DOUBLE, destRank, 101, MPI_COMM_WORLD);
          MPI_Send(&toSendPos[head_offset], numGalToRanks[destRank], MPI_INT, destRank, 102, MPI_COMM_WORLD);
        }
      }
      
      if((thisRank == destRank) && (head == destRank))
      {
        /* (send) copy stuff within rank */
        recvNumGalOnRank[head] = numGalToRanks[destRank];
        
        counterRecvNumGal += recvNumGalOnRank[head];
        (*toRecvNion) = realloc((*toRecvNion), counterRecvNumGal * sizeof(double));
        (*toRecvPos) = realloc((*toRecvPos), counterRecvNumGal * sizeof(int32_t));
        
        for(int i=0; i<recvNumGalOnRank[head]; i++)
        {
            (*toRecvNion)[counterRecvNumGal - recvNumGalOnRank[head] + i] = toSendNion[head_offset + i];
            (*toRecvPos)[counterRecvNumGal - recvNumGalOnRank[head] + i] = toSendPos[head_offset + i];
        }
      }
      
      head_offset += numGalToRanks[destRank];
    }
    
    if(thisRank != head)
    {
      /* receive stuff from head */
      MPI_Recv(&recvNumGalOnRank[head], 1, MPI_INT, head, 100, MPI_COMM_WORLD, &status);
      
      counterRecvNumGal += recvNumGalOnRank[head];
      if(recvNumGalOnRank[head] > 0)
      {
        (*toRecvNion) = realloc((*toRecvNion), counterRecvNumGal * sizeof(double));
        (*toRecvPos) = realloc((*toRecvPos), counterRecvNumGal * sizeof(int32_t));
              
        MPI_Recv(&(*toRecvNion)[counterRecvNumGal - recvNumGalOnRank[head]], recvNumGalOnRank[head], MPI_DOUBLE, head, 101, MPI_COMM_WORLD, &status);
        MPI_Recv(&(*toRecvPos)[counterRecvNumGal - recvNumGalOnRank[head]], recvNumGalOnRank[head], MPI_INT, head, 102, MPI_COMM_WORLD, &status);
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  numGalOnThisRank = &recvNumGalOnRank;
}
#endif

/*=======================================================================================*/

/*-----------------------------------------------------*/
/* COPYING GALAXY PROPERTIES TO NION LIST TO BE MAPPED */
/*-----------------------------------------------------*/

void copy_gal_to_nion(gal_t *thisGal, dconfObj_t simParam, domain_t *thisDomain, double factor, nion_t *thisNion)
{
  int32_t nbins = thisDomain->nbins;
  double inv_boxsize = simParam->inv_boxsize;
  
  if(thisNion->numGalWritten == thisNion->numGal)
  {
    reallocNion(&thisNion, (int32_t)(factor * thisNion->numGal));
  }
  
  int32_t index = thisNion->numGalWritten;
  int32_t x = 0, y = 0, z = 0;
  thisNion->nion[index] = 1.e-50 * get_nion_for_model(thisGal, simParam);//thisGal->Mvir * 1.e45;
  x = thisGal->pos[0] * inv_boxsize * nbins;  // assumes pos[i] to arange between 0 and 1
  y = thisGal->pos[1] * inv_boxsize * nbins;
  z = thisGal->pos[2] * inv_boxsize * nbins;
  if(x >= nbins)
    x = nbins - 1;
  if(y >= nbins)
    y = nbins - 1;
  if(z >= nbins)
    z = nbins - 1;
  thisNion->pos[index] = z*nbins*nbins + y*nbins + x;
  thisNion->rank[index] = calc_rank_from_pos(thisDomain, thisGal->pos[2]*inv_boxsize);
  
  thisNion->numGalWritten++;
}
