#ifndef TEST_H
#define TEST_H

#include "gal_gtree.h"

double cutandresort(gtree_t **theseGtrees, int numGtrees, int outSnap, gtree_t ***newGtrees);
int getNumberOfRoot(gtree_t *thisGtree, int outSnap);
void add_index(int ***listEquals, int *sizeListEquals, int index);
int search_index(int index, int **listEquals);
int isNotInThisList(int index, int *thisList);
void mergeLists(int ***listEquals, int *sizeListEquals, int listOne, int listTwo);
void addGalToTree(gtree_t ***thisGtree, int position, gal_t thisGal);
void addTrees(gtree_t ***theseNewGtrees, int *numberNewGtrees, gtree_t **new_trees, int sizeListEquals);

#endif 
