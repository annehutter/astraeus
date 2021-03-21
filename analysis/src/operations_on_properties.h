#ifndef OPERATIONS_ON_PROPERTIES_H
#define OPERATIONS_ON_PROPERTIES_H

int check_operation_present(char *propertyName);
void get_delimiters_in_propertyName(char *propertyName, int *numDelim, char **listDelim);
void extract_properties(char *propertyName, int *numProperties, char ***listPropertyNames);
void multiply_properties(int num, double *property1, double* property2);
void divide_properties(int num, double *property1, double* property2);

int check_MvirThreshold_present(char *propertyName);
char *get_propertyName(char *propertyName);
int get_num_reducedProperty(int numGal, int numSnaps, double *Mvir, double MvirThreshold);
double *reduceProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times, double *thisPropertyHistory, int *numGalSelectedProperty);
double *reducePropertyHistory(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times, double *thisProperty, int *numGalSelectedProperty);

double *selectPropertyWithMapping(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, double propertyLowLimit, double propertyUpLimit, int currSnap, float *times);
double *selectProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times);

double *getThisPropertyWithMapping(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, double propertyLowLimit, double propertyUpLimit, int currSnap, float *times, int *numGalProperty);
double *getThisProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times, int *numGalProperty);
double *getThisPropertyHistory(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times, int *numGalProperty);

#endif