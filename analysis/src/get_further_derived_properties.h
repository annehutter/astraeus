#ifndef GET_FURTHER_DERIVED_PROPERTIES_H
#define GET_FURTHER_DERIVED_PROPERTIES_H

void prepare_stellarmasshistory(double *stellarmasshistory, int numGal, int currSnap);

int check_derived_property_present(char *propertyName);

double get_max_of_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times);

double *get_derived_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times);

double *get_ionEmissivity_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times);
double *get_SNejection_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times);

double *get_fesc_SNejection_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times);
double *get_fesc_property(dconfObj_t simParam, int numTrees, int *sizeListEquals);

#endif