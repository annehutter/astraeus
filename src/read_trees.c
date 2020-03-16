#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"

#include "halo_tree.h"
#include "read_trees.h"

tree_t *read_tree(FILE *f)
{
  tree_t *newTree = NULL;
  int32_t numHalos = 0;
  
  /* READ HOW MANY HALOS ARE IN THIS TREE */
  fread(&numHalos, sizeof(int32_t), 1 ,f);
    
  newTree = initTree(numHalos);
  newTree->numHalos = numHalos;
  
  fread(newTree->halos, sizeof(halo_t), numHalos, f);
  
  return newTree;
}

int32_t read_trees_in_file(char *fileName, tree_t ***thisTreeList, int offset)
{
  tree_t **theseTrees = *thisTreeList;
  int32_t numTrees = 0;
  int32_t numTreesTmp = 0;
  
  FILE *f = NULL;
  
  if( (f=fopen(fileName, "r")) == NULL)
  {
    fprintf(stderr, "Could not open file\n");
    exit(EXIT_FAILURE);
  }

  /* READING TREES */
  
//   while(!feof(f))
//   {
    fread(&numTreesTmp, sizeof(int32_t), 1, f);
     
    printf("numTrees read = %d\n", numTreesTmp);
    theseTrees = realloc(theseTrees, sizeof(tree_t) * (numTrees + numTreesTmp + offset));
    for(int tree=0; tree<numTreesTmp; tree++)
    {
      theseTrees[offset + numTrees + tree] = read_tree(f);
    }
    numTrees += numTreesTmp;
//   }

  *thisTreeList = theseTrees;
  
  fclose(f);
  
  return numTrees;
}

int32_t *get_num_files_to_read(int numFiles, int thisRank, int size)
{  
  int32_t *numFilesToRead = allocate_array_int32_t(size, "numFilesToRead");

  /* DETERMINE HOW MANY FILES EACH PROCESSOR NEEDS TO READ */
#ifdef MPI
  int32_t numFilesOnMostRanks = numFiles / size;
  int32_t numFilesThisRank = numFilesOnMostRanks;
  
  int32_t numFilesAdd = numFiles - numFilesOnMostRanks * size;
  if(thisRank < numFilesAdd)
    numFilesThisRank = numFilesOnMostRanks + 1;
  
  MPI_Allgather(&numFilesThisRank, 1, MPI_INT, numFilesToRead, 1, MPI_INT, MPI_COMM_WORLD);
#else
  if(size == 1) numFilesToRead[0] = numFiles;
  else printf("Running in serial but size != 1. Something is seriously wrong!\n");
#endif
  
  return numFilesToRead;
}

char *get_file_name(char *baseName, int32_t file, int32_t *numFilesToRead, int thisRank)
{
  char *thisFileName = NULL;
  char fileEnding[10];
  int32_t offset = 0;
  
  for(int rank=0; rank<thisRank; rank++) 
    offset += numFilesToRead[rank];
  offset += file;

  sprintf(fileEnding, "_%d.dat", offset);
  thisFileName = concat(baseName, fileEnding);
  
  return thisFileName;
}

/* reading in routine that can deal with reading multiple files (each different on eavch processor) */
int32_t read_trees_in_all_files(char *baseName, int numFiles, tree_t ***thisTreeList, int thisRank, int size)
{
  tree_t **theseTrees = *thisTreeList;
  int32_t numTrees = 0;
  int32_t offset = 0;
  char *fileName = NULL;
  
  int32_t *numFilesToRead = get_num_files_to_read(numFiles, thisRank, size);

  for(int file=0; file<numFilesToRead[thisRank]; file++)
  {
    fileName = get_file_name(baseName, file, numFilesToRead, thisRank);
    printf("rank %d: reading file '%s'\n", thisRank, fileName);
    offset = numTrees;
    numTrees += read_trees_in_file(fileName, &theseTrees, offset);
    free(fileName);
  }
  free(numFilesToRead);
  
  *thisTreeList = theseTrees;
  
  return numTrees;
}
