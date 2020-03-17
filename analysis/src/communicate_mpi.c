#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "sort.h"
#include "get_grid_properties.h"
#include "communicate_mpi.h"

/*----------------------------------------------------------------*/
/* CREATING & SORTING ARRAY FOR SENDING */
/*----------------------------------------------------------------*/

int *create_rankArray_3Dgrid(int size, int numEntries, float *posIndex, int gridsize)
{
  int *rankArray = allocate_array_int(numEntries, "rankArray");
  
  float factor = (float)size / (float)(gridsize*gridsize*gridsize);
  for(int i=0; i<numEntries; i++)
  {
    rankArray[i] = ((int)(posIndex[i] * factor) < size ? (int)(posIndex[i] * factor) : size-1);
  }
  
  return rankArray;
}

/* Creating index array for sorting */
int *create_sorted_indexArray(int numEntries, int *arrayToSort)
{
  int32_t *indexArray = allocate_array_int32_t(numEntries, "indexArray");
  for(int i=0; i<numEntries; i++)
    indexArray[i] = i;
  
  quicksort(arrayToSort, 0, numEntries-1, indexArray);
  
  return indexArray;
}

int *create_sortedArray_int(int numEntries, int *thisArray, int *rankArray, int **indexArray)
{
  *indexArray = create_sorted_indexArray(numEntries, rankArray);
  int *sortedArray = allocate_array_int(numEntries, "sortedArray");
  
  for(int i=0; i<numEntries; i++)
  {
    sortedArray[i] = thisArray[(*indexArray)[i]];
  }
  
  return sortedArray;
}

float *create_sortedArray_float(int numEntries, float *thisArray, int *rankArray, int **indexArray)
{
  *indexArray = create_sorted_indexArray(numEntries, rankArray);
  float *sortedArray = allocate_array_float(numEntries, "sortedArray");
  
  for(int i=0; i<numEntries; i++)
  {
    sortedArray[i] = thisArray[(*indexArray)[i]];
  }
  
  return sortedArray;
}

double *create_sortedArray_double(int numEntries, double *thisArray, int *rankArray, int **indexArray)
{
  *indexArray = create_sorted_indexArray(numEntries, rankArray);
  double *sortedArray = allocate_array_double(numEntries, "sortedArray");
  
  for(int i=0; i<numEntries; i++)
  {
    sortedArray[i] = thisArray[(*indexArray)[i]];
  }
  
  return sortedArray;
}

int32_t *create_numToSendToRanksArray(int size, int numEntries, int *rankArray, int *indexArray)
{
  int32_t *numToSendToRanks = allocate_array_int32_t(size, "numToSendToRanks");
  
  int counter = 0;
  int rank = 0;
  for(int i=0; i<numEntries; i++)
  {
//     if(i<3 || i>numEntries-3) printf("rankArray[%d] = %d\n", i, rankArray[i]);
    if(rankArray[i] == rank)
    {
      counter++;
    }
    else
    {
      numToSendToRanks[rank] = counter;
      counter = 1;
      rank++;
    }
  }
  numToSendToRanks[rank] = counter;
  
  int sum = 0;
  for(int i=0; i<size; i++)
    sum += numToSendToRanks[i];
  
  if(sum != numEntries)
  {
    for(int i=0; i<size; i++)
      printf("num[%d] = %d\t", i, numToSendToRanks[i]);
    printf("sum = %d\t numEntries = %d\n", sum, numEntries);
  }
  assert(sum == numEntries);
  
  return numToSendToRanks;
}

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/

int get_numEntries_numOnThisRank(int size, int32_t *numOnThisRank)
{
  int sum = 0;
  
  for(int i=0; i<size; i++)
    sum += numOnThisRank[i];
  
  return sum;
}

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/

float *create_resortedArray(int numEntries, float *sortedArray, int *indexArray)
{
  float *resortedArray = allocate_array_float(numEntries, "resortedArray");
  
  for(int i=0; i<numEntries; i++)
  {
    resortedArray[indexArray[i]] = sortedArray[i];
  }
  
  return resortedArray;
}

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/

#ifdef MPI
float *communicate_array_float(int size, int thisRank, int numEntries, float *posIndex, int gridsize, float *thisArray, float *thisGrid)
{
  /*--------------------------------*/
  /* prepare all necessary arrays */
  /*--------------------------------*/

  int *indexArray = NULL;
  
  /* create rankArray */
  int *rankArray = create_rankArray_3Dgrid(size, numEntries, posIndex, gridsize);
  
  /* create sortedArray = thisArray sorted by ranks */
  float *sortedArray = create_sortedArray_float(numEntries, thisArray, rankArray, &indexArray);
  
  /* create array with number of galaxies to send to ranks */
  int32_t *numOnRanks = create_numToSendToRanksArray(size, numEntries, rankArray, indexArray);
  
  /* arrays for receiving */
  int32_t *recvNumOnRanks = allocate_array_int32_t(size, "recvNumOnRanks");
  int32_t numRecvEntries = 0;
  float *recvArray = NULL;
  float *recvNewArray = NULL;
  float *resortedArray = NULL;
    
  /*--------------------------------*/
  /* sending forward */
  /*--------------------------------*/  
  
  send_recv_array_float(size, thisRank, numOnRanks, &recvNumOnRanks, sortedArray, &recvArray);  
  free(sortedArray);
  sortedArray = NULL;
  
  /*--------------------------------*/
  /* create new array */
  /*--------------------------------*/  
  
  numRecvEntries = get_numEntries_numOnThisRank(size, recvNumOnRanks);
  recvNewArray = allocate_array_float(numRecvEntries, "recvNewArray");
  get_grid_values(size, thisRank, numRecvEntries, &recvNewArray, recvArray, gridsize, thisGrid);
  free(recvArray);
  
  /*--------------------------------*/
  /* sending backward */
  /*--------------------------------*/  
  
  send_recv_array_float(size, thisRank, recvNumOnRanks, &numOnRanks, recvNewArray, &sortedArray);
  free(recvNewArray);

  /*--------------------------------*/
  /* resorting */
  /*--------------------------------*/ 
  
  resortedArray = create_resortedArray(numEntries, sortedArray, indexArray);
  
  /*--------------------------------*/
  /* deallocation */
  /*--------------------------------*/
  free(indexArray);
  free(rankArray);
  free(sortedArray);
  free(numOnRanks);
  free(recvNumOnRanks);
  
  return resortedArray;
}
#endif

/*----------------------------------------------------------------*/
/* SEND & RECEIVE SORTED ARRAY ACROSS PROCESSORS */
/*----------------------------------------------------------------*/

#ifdef MPI
void send_recv_array_double(int size, int thisRank, int32_t *numToSendToRanks, int32_t **numOnThisRank, double *toSendArray, double **toRecvArray)
/* numToSendToRanks: array[size] that holds the number of entries to send to each rank 
   numOnThisRank: *array[size] that holds number which have been received from each rank 
   toSendArray: rank ordered array to be sent
   toRecvArray: pointer to array received from all ranks */
{
  MPI_Status status;
  
  int32_t head_offset = 0;
  
  int32_t *recvNumOnRank = *numOnThisRank;
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
        MPI_Send(&numToSendToRanks[destRank], 1, MPI_INT, destRank, 400, MPI_COMM_WORLD);
        
        if(numToSendToRanks[destRank] > 0)
        {            
          MPI_Send(&toSendArray[head_offset], numToSendToRanks[destRank], MPI_DOUBLE, destRank, 401, MPI_COMM_WORLD);
        }
      }
      
      if((thisRank == destRank) && (head == destRank))
      {
        /* (send) copy stuff within rank */
        recvNumOnRank[head] = numToSendToRanks[destRank];
        
        counterRecvNumGal += recvNumOnRank[head];
        (*toRecvArray) = realloc((*toRecvArray), counterRecvNumGal * sizeof(double));
        
        for(int i=0; i<recvNumOnRank[head]; i++)
        {
            (*toRecvArray)[counterRecvNumGal - recvNumOnRank[head] + i] = toSendArray[head_offset + i];
        }
      }
      
      head_offset += numToSendToRanks[destRank];
    }
    
    if(thisRank != head)
    {
      /* receive stuff from head */
      MPI_Recv(&recvNumOnRank[head], 1, MPI_INT, head, 400, MPI_COMM_WORLD, &status);
      
      counterRecvNumGal += recvNumOnRank[head];
      if(recvNumOnRank[head] > 0)
      {
        (*toRecvArray) = realloc((*toRecvArray), counterRecvNumGal * sizeof(double));
              
        MPI_Recv(&(*toRecvArray)[counterRecvNumGal - recvNumOnRank[head]], recvNumOnRank[head], MPI_DOUBLE, head, 401, MPI_COMM_WORLD, &status);
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  numOnThisRank = &recvNumOnRank;
}
#endif

#ifdef MPI
void send_recv_array_float(int size, int thisRank, int32_t *numToSendToRanks, int32_t **numOnThisRank, float *toSendArray, float **toRecvArray)
/* numToSendToRanks: array[size] that holds the number of entries to send to each rank 
   numOnThisRank: *array[size] that holds number which have been received from each rank 
   toSendArray: rank ordered array to be sent
   toRecvArray: pointer to array received from all ranks */
{
  MPI_Status status;
  
  int32_t head_offset = 0;
  
  int32_t *recvNumOnRank = *numOnThisRank;
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
        MPI_Send(&numToSendToRanks[destRank], 1, MPI_INT, destRank, 200, MPI_COMM_WORLD);
        
        if(numToSendToRanks[destRank] > 0)
        {            
          MPI_Send(&toSendArray[head_offset], numToSendToRanks[destRank], MPI_FLOAT, destRank, 201, MPI_COMM_WORLD);
        }
      }
      
      if((thisRank == destRank) && (head == destRank))
      {
        /* (send) copy stuff within rank */
        recvNumOnRank[head] = numToSendToRanks[destRank];
        
        counterRecvNumGal += recvNumOnRank[head];
        (*toRecvArray) = realloc((*toRecvArray), counterRecvNumGal * sizeof(float));
        
        for(int i=0; i<recvNumOnRank[head]; i++)
        {
            (*toRecvArray)[counterRecvNumGal - recvNumOnRank[head] + i] = toSendArray[head_offset + i];
        }
      }
      
      head_offset += numToSendToRanks[destRank];
    }
    
    if(thisRank != head)
    {
      /* receive stuff from head */
      MPI_Recv(&recvNumOnRank[head], 1, MPI_INT, head, 200, MPI_COMM_WORLD, &status);
      
      counterRecvNumGal += recvNumOnRank[head];
      if(recvNumOnRank[head] > 0)
      {
        (*toRecvArray) = realloc((*toRecvArray), counterRecvNumGal * sizeof(float));
              
        MPI_Recv(&(*toRecvArray)[counterRecvNumGal - recvNumOnRank[head]], recvNumOnRank[head], MPI_FLOAT, head, 201, MPI_COMM_WORLD, &status);
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  numOnThisRank = &recvNumOnRank;
}
#endif

#ifdef MPI
void send_recv_array_int(int size, int thisRank, int32_t *numToSendToRanks, int32_t **numOnThisRank, int *toSendArray, int **toRecvArray)
/* numToSendToRanks: array[size] that holds the number of entries to send to each rank 
   numOnThisRank: *array[size] that holds number which have been received from each rank 
   toSendArray: rank ordered array to be sent
   toRecvArray: pointer to array received from all ranks */
{
  MPI_Status status;
  
  int32_t head_offset = 0;
  
  int32_t *recvNumOnRank = *numOnThisRank;
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
        MPI_Send(&numToSendToRanks[destRank], 1, MPI_INT, destRank, 100, MPI_COMM_WORLD);
        
        if(numToSendToRanks[destRank] > 0)
        {            
          MPI_Send(&toSendArray[head_offset], numToSendToRanks[destRank], MPI_INT, destRank, 101, MPI_COMM_WORLD);
        }
      }
      
      if((thisRank == destRank) && (head == destRank))
      {
        /* (send) copy stuff within rank */
        recvNumOnRank[head] = numToSendToRanks[destRank];
        
        counterRecvNumGal += recvNumOnRank[head];
        (*toRecvArray) = realloc((*toRecvArray), counterRecvNumGal * sizeof(int));
        
        for(int i=0; i<recvNumOnRank[head]; i++)
        {
            (*toRecvArray)[counterRecvNumGal - recvNumOnRank[head] + i] = toSendArray[head_offset + i];
        }
      }
      
      head_offset += numToSendToRanks[destRank];
    }
    
    if(thisRank != head)
    {
      /* receive stuff from head */
      MPI_Recv(&recvNumOnRank[head], 1, MPI_INT, head, 100, MPI_COMM_WORLD, &status);
      
      counterRecvNumGal += recvNumOnRank[head];
      if(recvNumOnRank[head] > 0)
      {
        (*toRecvArray) = realloc((*toRecvArray), counterRecvNumGal * sizeof(int));
              
        MPI_Recv(&(*toRecvArray)[counterRecvNumGal - recvNumOnRank[head]], recvNumOnRank[head], MPI_INT, head, 101, MPI_COMM_WORLD, &status);
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  numOnThisRank = &recvNumOnRank;
}
#endif