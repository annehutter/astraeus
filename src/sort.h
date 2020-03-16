#ifndef SORT_H
#define SORT_H

void quicksort(int *array, int low, int high, int *indexArray);
int partition(int *array, int low, int high, int *indexArray);
void swap(int *array, int i, int j);

#endif