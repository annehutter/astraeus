#ifndef READ_H
#define READ_H

FILE *openToRead(char *infilename);

outgtree *read_tree(FILE *f);
int read_trees_in_file(char *fileName, outgtree ***thisTreeList);

#endif
