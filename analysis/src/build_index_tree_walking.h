#ifndef BUILD_INDEX_TREE_WALKING_H
#define BUILD_INDEX_TREE_WALKING_H

void add_index(int ***listEquals, int *sizeListEquals, int index);
int search_index(int index, int **listEquals);
int is_not_in_thisList(int index, int *thisList);

int *init_toBeMerged();
void write_in_toBeMergedList(int indexGal, int indexDescGal, int **toBeMerged);
int is_index_in_toBeMergedList(int indexDescGal, int **toBeMerged);
int get_index_toBeMergedList(int indexDescGal, int **toBeMerged);
void print_toBeMergedList(int *toBeMerged);

void mergeLists(int ***listEquals, int *sizeListEquals, int listOne, int listTwo);

void printlistequals(int **listEquals, int numberlist);

void get_listEquals(outgtree_t **theseGtrees, int numGtrees, int **index, int ***listEquals, int *sizeListEquals, int **listMerged, int prevSnap, int currSnap);

#endif