#ifndef TREE_WALKING
#define TREE_WALKING

int getmax(int *listOfInt, int length);
void printlistequals(int **listEquals, int numberlist);
void getListEquals(outgtree_t **theseGtrees, int numGtrees, int **index, int ***listEquals, int *sizeListEquals, int prevSnap, int currentSnap);
void add_index(int ***listEquals, int *sizeListEquals, int index);
int search_index(int index, int **listEquals);
int isNotInThisList(int index, int *thisList);
void mergeLists(int ***listEquals, int *sizeListEquals, int listOne, int listTwo);

#endif