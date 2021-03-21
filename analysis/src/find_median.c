#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "find_median.h"

void swap_float(float *a, float *b)
{
  float tmp = *a;
  *a = *b;
  *b = tmp;
}

// Have an array of length num on each processor & each array needs to be sorted according to a pivot value
long int sortArray_with_pivot(float *sortArray, long int num, float pivot, long int low, long int high)
{
  assert(high <= num);
  long int i = (low - 1);
  
  for(long int j=low; j<high; j++)
  {
    if(sortArray[j] < pivot)
    {
      i++;
      swap_float(&(sortArray[i]), &(sortArray[j]));
    }
  }
//   swap_float(&(sortArray[i+1]), &(sortArray[pivot]));
  return i+1;
}

float get_median(float *sortArray, long int numArray, int thisRank, int size, float minPivot, float maxPivot)
{
  /* LOCAL ON EACH PROCESSOR */
  long int cut = 0;
  float pivot = 0., newPivot = 0.;
  long int low = 0;
  long int high = numArray;
  
  long int numSum = 0;
  long int cutSum = 0;
  long int *lowOnRanks = NULL;
  long int *highOnRanks = NULL;
  long int *cutOnRanks = NULL;
  long int *numOnRanks = allocate_array_long_int(size, "numOnRanks");

  /* ONLY ON MAIN PROCESSOR */
  if(thisRank == 0)
  {
    lowOnRanks = allocate_array_long_int(size, "lowOnRanks");
    highOnRanks = allocate_array_long_int(size, "highOnRanks");
    cutOnRanks = allocate_array_long_int(size, "cutOnRanks");
  }
  
  /* communicate array sizes to main processor */
#ifdef MPI
  MPI_Allgather(&numArray, 1, MPI_LONG, numOnRanks, 1, MPI_LONG, MPI_COMM_WORLD);
#else
  numArray = numOnRanks[0];
#endif  
  for(int rank=0; rank<size; rank++)
  {
    numSum += numOnRanks[rank];
    if(thisRank == 0) highOnRanks[rank] = numOnRanks[rank];
  }
  
  /* ----- LOOP ------ */
  
  if(numSum <= 2 && numSum > 0)
  {
    for(int rank=0; rank<size; rank++)
    {
      if(numOnRanks[rank] > 0)
      {
#ifdef MPI
        if(thisRank == rank && rank != 0)
        {
          MPI_Ssend(sortArray, numOnRanks[rank], MPI_FLOAT, 0, 100, MPI_COMM_WORLD);
        }
#endif
        if(thisRank == 0)
        {
          if(rank != 0)
          {
            float *recvArray = allocate_array_float(numOnRanks[rank], "recvArray");
#ifdef MPI
            MPI_Recv(recvArray, numOnRanks[rank], MPI_FLOAT, rank, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif
            for(int i=0; i<numOnRanks[rank]; i++)
            {
                pivot += recvArray[i];
            }
            free(recvArray);
          }
          else
          {
            for(int i=0; i<numOnRanks[rank]; i++)
            {
                pivot += sortArray[i];
            }
          }
        }
      }
    }
    if(numSum > 0)
      pivot = pivot / (float)numSum;
  }
  else
  {
    int loop = 0;
    while(labs(numSum/2 - cutSum) > 0.5)
    {
      // choose & communicate pivot
      if(thisRank == 0)
      {
        float tmp = ((float)rand() / (float)RAND_MAX);
        newPivot = minPivot + tmp * (maxPivot - minPivot);
      }
      
#ifdef MPI
      MPI_Bcast(&pivot, 1, MPI_LONG, 0, MPI_COMM_WORLD);
      MPI_Bcast(&newPivot, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif

      if(loop ==0 || (newPivot/pivot < 0.999 || newPivot/pivot > 1.001))
        pivot = newPivot;
      else
      {
        break;
      }
      
      // sort array according to pivot (local)
      cut = sortArray_with_pivot(sortArray, numArray, pivot, low, high);

      // cut provides index where higher half starts = number of lower half
#ifdef MPI
      MPI_Gather(&cut, 1, MPI_LONG, cutOnRanks, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#else
      cut = cutOnRanks[0];
#endif

      // on main processor
      if(thisRank == 0)
      {
        cutSum = 0;
        for(int rank=0; rank<size; rank++)
        {
          cutSum += cutOnRanks[rank];
        }
              
        for(int rank=0; rank<size; rank++)
        {
          if(cutSum < numSum/2)
          {
            /* disregard low values */
            lowOnRanks[rank] = cutOnRanks[rank];
            minPivot = pivot;
          }
          else if(cutSum > numSum/2)
          {
            /* disregard high values */
            highOnRanks[rank] = cutOnRanks[rank];
            maxPivot = pivot;
          }
        }
      }
    
#ifdef MPI
      MPI_Scatter(lowOnRanks, 1, MPI_LONG, &low, 1, MPI_LONG, 0, MPI_COMM_WORLD);
      MPI_Scatter(highOnRanks, 1, MPI_LONG, &high, 1, MPI_LONG, 0, MPI_COMM_WORLD);    
      MPI_Bcast(&cutSum, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#else
      low = lowOnRanks[0];
      high = highOnRanks[0];
#endif

      loop++;
    }
  }
  
#ifdef MPI
  MPI_Bcast(&pivot, 1, MPI_LONG, 0, MPI_COMM_WORLD);
#endif
  
  if(numOnRanks != NULL) free(numOnRanks);
  if(cutOnRanks != NULL) free(cutOnRanks);
  if(lowOnRanks != NULL) free(lowOnRanks);
  if(highOnRanks != NULL) free(highOnRanks);
  
  return pivot;
}