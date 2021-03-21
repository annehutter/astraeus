#ifndef READ_TREES_H
#define READ_TREES_H

void print_tree(tree_t *thisTree);

tree_t *read_tree(FILE *f);
int32_t read_trees_in_file(char *fileName, tree_t ***thisTreeList, int offset);

int32_t *get_num_files_to_read(int numFiles, int thisRank, int size);
char *get_file_name(char *baseName, int32_t file, int32_t *numFilesToRead, int thisRank);
int32_t read_trees_in_all_files(char *baseName, int numFiles, tree_t ***thisTreeList, int thisRank, int size);

#endif