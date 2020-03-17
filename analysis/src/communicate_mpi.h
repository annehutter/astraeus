#ifndef COMMUNICATE_MPI_H
#define COMMUNICATE_MPI_H

int *create_rankArray_3Dgrid(int size, int numEntries, float *posIndex, int gridsize);
int *create_sorted_indexArray(int numEntries, int *arrayToSort);
int *create_sortedArray_int(int numEntries, int *thisArray, int *rankArray, int **indexArray);
float *create_sortedArray_float(int numEntries, float *thisArray, int *rankArray, int **indexArray);
double *create_sortedArray_double(int numEntries, double *thisArray, int *rankArray, int **indexArray);
int32_t *create_numToSendToRanksArray(int size, int numEntries, int *rankArray, int *indexArray);

int get_numEntries_numOnThisRank(int size, int32_t *numOnThisRank);
float *create_resortedArray(int numEntries, float *sortedArray, int *indexArray);
#ifdef MPI
float *communicate_array_float(int size, int thisRank, int numEntries, float *posIndex, int gridsize, float *thisArray,  float *thisGrid);
#endif

#ifdef MPI
void send_recv_array_double(int size, int thisRank, int32_t *numToSendToRanks, int32_t **numOnThisRank, double *toSendArray, double **toRecvArray);
void send_recv_array_float(int size, int thisRank, int32_t *numToSendToRanks, int32_t **numOnThisRank, float *toSendArray, float **toRecvArray);
void send_recv_array_int(int size, int thisRank, int32_t *numToSendToRanks, int32_t **numOnThisRank, int *toSendArray, int **toRecvArray);
#endif

#endif