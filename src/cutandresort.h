#ifndef TEST_H
#define TEST_H

#include "gal_gtree.h"

double cutandresort(gtree_t **theseGtrees, int numGtrees, int outSnap, gtree_t ***newGtrees);
int search_index(int *arrayToSearch, int value);
void reattribute_localID(gtree_t **thisGtrees, int numberTrees);
int get_number_of_root(gtree_t *thisGtree, int outSnap);
void add_index(int ***listEquals, int *sizeListEquals, int index);
int search_branch(int index, int **listEquals);
int is_not_in_this_list(int index, int *thisList);
void merge_lists(int ***listEquals, int *sizeListEquals, int listOne, int listTwo);
void add_gal_to_tree(gtree_t ***thisGtree, int position, gal_t thisGal);
void add_trees(gtree_t ***theseNewGtrees, int *numberNewGtrees, gtree_t **new_trees, int sizeListEquals);

#endif 
