#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "sort.h"

/* QUICKSORT ALGORITHM */

void quicksort(int *array, int low, int high, int *indexArray)
{
  int p = 0;
  
  if(low < high)
  {
    p = partition(array, low, high, indexArray);
    quicksort(array, low, p-1, indexArray);
    quicksort(array, p+1, high, indexArray);
  }
}

int partition(int *array, int low, int high, int *indexArray)
{
  int pivot = array[high];
  int i = low;
  
  for(int j=low; j<=high-1; j++)
  {
    if(array[j] < pivot)
    {
      swap(array, i, j);
      swap(indexArray, i, j);
      i++;
    }
  }
  swap(array, i, high);
  swap(indexArray, i, high);
  
  return i;
}

void swap(int *array, int i, int j)
{
  int tmpArray = array[i];
  array[i] = array[j];
  array[j] = tmpArray;
}