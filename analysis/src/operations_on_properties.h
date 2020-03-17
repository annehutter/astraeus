#ifndef OPERATIONS_ON_PROPERTIES_H
#define OPERATIONS_ON_PROPERTIES_H

int check_operation_present(char *propertyName);
void get_delimiters_in_propertyName(char *propertyName, int *numDelim, char **listDelim);
void extract_properties(char *propertyName, int *numProperties, char ***listPropertyNames);
void multiply_properties(int num, float *property1, float* property2);
void divide_properties(int num, float *property1, float* property2);
float *getGridProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times);
float *getThisProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times);
float *getThisPropertyHistory(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times);

#endif