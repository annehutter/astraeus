#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "utils.h"
#include "outgal.h"
#include "read.h"

FILE *openToRead(char *infilename)
{
    FILE *infile = NULL;
    if((infile = fopen(infilename, "r")) == NULL)
    {
        printf("Can't open file at location %s\n", infilename);
        exit(EXIT_FAILURE);
    }
    return infile;
}

outgtree *read_tree(FILE *f)
{
    outgtree *newTree = NULL;
    int32_t numGals = 0;
    
    /* READ HOW MANY HALOS ARE IN THIS TREE */
    fread(&numGals, sizeof(int32_t), 1 ,f);
    newTree = initOutGtree(numGals);
    newTree->numGal = numGals;
    
    fread(newTree->outgalaxies, sizeof(outgal_t), numGals, f);
    
    return newTree;
}


int read_trees_in_file(char *fileName, outgtree ***thisTreeList)
{
    outgtree **theseTrees = *thisTreeList;
    int32_t numTrees = 0;
    FILE *f = NULL;
    
    f = openToRead(fileName);

    /* READING TREES */
    fread(&numTrees, sizeof(int32_t), 1, f);
    printf("numTrees read = %d\n", numTrees);
    
    theseTrees = malloc(sizeof(outgtree*)* numTrees);
    for(int tree=0; tree<numTrees; tree++)
            theseTrees[tree] = read_tree(f);
    *thisTreeList = theseTrees;

    fclose(f);
  
    return numTrees;
}
