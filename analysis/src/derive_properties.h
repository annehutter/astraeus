#ifndef DERIVE_PROPERTIES_H
#define DERIVE_PROPERTIES_H

int get_numGal_at_snap(int *sizeListEquals, int numTrees);
float *getProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times);
float *getPropertyHistory(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times);

#endif
