#ifndef FIND_MEDIAN_H
#define FIND_MEDIAN_H

void swap_float(float *a, float *b);
long int sortArray_with_pivot(float *sortArray, long int num, float pivot, long int low, long int high);
float get_median(float *sortArray, long int numArray, int thisRank, int size, float minPivot, float maxPivot);

#endif