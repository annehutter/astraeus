#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "dconfObj.h"
#include "outgal.h"
#include "statistics.h"

#include "build_index_tree_walking.h"

void printlistequals(int **listEquals, int numberlist)
{
        for(int i=0; i<numberlist; i++)
        {
                printf("%d entries: ", listEquals[i][0]);
                printf("[");
                for(int j=1; j<listEquals[i][0]+1; j++)
                        printf(" %d ",listEquals[i][j]);
                printf("]  ");
        }
        printf("\n");
}

void add_index(int ***listEquals, int *sizeListEquals, int index)
{
  /* increasing the number of entries by one */
  (*sizeListEquals)++;
  
  /* reallocating listEquals[tree] */
  int **ptr_tmp = realloc(*listEquals, (*sizeListEquals)*sizeof(int*));
  if(ptr_tmp == NULL)
  {
    printf("realloc failed, exiting !\n");
    free(*listEquals);
    exit(0);
  }
  else
    *listEquals = ptr_tmp;
  
  (*listEquals)[*sizeListEquals-1] = malloc(2*sizeof(int));
  (*listEquals)[*sizeListEquals-1][0] = 1; //first element of the list = size of the list - 1
  (*listEquals)[*sizeListEquals-1][1] = index;
}

int search_index(int index, int **listEquals)
{
  int wantedList = 0;

  while(is_not_in_thisList(index, listEquals[wantedList]))
  {
//     printf("listEquals[%d][0] = %d\twantedList = %d\n", wantedList, listEquals[wantedList][0], wantedList);
    wantedList++;
  }
  
  return wantedList;
}

int is_not_in_thisList(int index, int *thisList)
{
  int isNot = 1;
  for(int i=1; i<thisList[0]+1; i++)
    if(index == thisList[i])
      isNot = 0;
  
  return isNot;
}

int *init_toBeMerged()
{
  int *toBeMerged = allocate_array_int(1, "toBeMerged");
  toBeMerged[0] = 0;
  
  return toBeMerged;
}

void write_in_toBeMergedList(int indexGal, int indexDescGal, int **toBeMerged)
{
  int numEntries = (*toBeMerged)[0] + 1;
  int *tmp = realloc(*toBeMerged, (numEntries*2 + 1) * sizeof(int));
  if(tmp == NULL)
  {
    fprintf(stderr, "reallocating of toBeMerged failed\n");
    exit(EXIT_FAILURE);
  }
  else
  *toBeMerged = tmp;
    
  (*toBeMerged)[0] = numEntries;
  (*toBeMerged)[2*numEntries-1] = indexGal;
  (*toBeMerged)[2*numEntries] = indexDescGal;
}

int is_index_in_toBeMergedList(int indexDescGal, int **toBeMerged)
{
  int isInList = 0;
  int numEntries = (*toBeMerged)[0];
  for(int i=0; i<numEntries; i++)
  {
    if((*toBeMerged)[2 + 2*i] == indexDescGal)
    {
      isInList = 1;
      break;
    }
  }
  
  return isInList;
}

int get_index_toBeMergedList(int indexDescGal, int **toBeMerged)
{
  int numEntries = (*toBeMerged)[0];
  int indexGal = -1;
  int entry = 0;
  for(int i=0; i<numEntries; i++)
  {
    if((*toBeMerged)[2 + 2*i] == indexDescGal)
    {
      indexGal = (*toBeMerged)[1 + 2*i];
      entry = i;
    }
  }
    
  numEntries--;
  (*toBeMerged)[0] = numEntries;
  
  for(int i=entry; i<numEntries; i++)
  {
    (*toBeMerged)[1 + 2*i] = (*toBeMerged)[1 + 2*(i+1)];
    (*toBeMerged)[2 + 2*i] = (*toBeMerged)[2 + 2*(i+1)];
  }
  
  int *tmp = realloc(*toBeMerged, (numEntries*2 + 1) * sizeof(int));
  if(tmp == NULL)
  {
    fprintf(stderr, "reallocating of toBeMerged failed\n");
    exit(EXIT_FAILURE);
  }
  else
    *toBeMerged = tmp;
  
  return indexGal;
}

void print_toBeMergedList(int *toBeMerged)
{
  int numEntries = toBeMerged[0];
  printf("numEntries = %d: ", numEntries);
  
  for(int i=0; i<numEntries; i++)
  {
    printf("%d -> %d,  ", toBeMerged[1 + 2*i], toBeMerged[2 + 2*i]);
  }
  printf("\n");
}

// merge two lists in the listEquals, free the second list and move the list with higher index accordingly
void merge_lists(int ***listEquals, int *sizeListEquals, int listOne, int listTwo)
{
  int **theseLists = *listEquals;
  int *ptr_tmp = NULL; 

  ptr_tmp = realloc(theseLists[listOne], sizeof(int)*(theseLists[listOne][0] + theseLists[listTwo][0]+1));
  
  if(ptr_tmp==NULL)
  {
    printf("realloc didn't work ! (test.c line 167)\n");
    exit(0);
  }
  else
    theseLists[listOne] = ptr_tmp;
          
  for(int i=1; i<theseLists[listTwo][0]+1; i++)
    theseLists[listOne][i+theseLists[listOne][0]] = theseLists[listTwo][i];
  
  theseLists[listOne][0] += theseLists[listTwo][0];       
  free(theseLists[listTwo]);
  for(int i=listTwo; i<(*sizeListEquals)-1; i++)
    theseLists[i] = theseLists[i+1];
          
  theseLists[(*sizeListEquals)-1] = NULL;
  
  (*sizeListEquals)--; 
}

void get_listEquals(outgtree_t **theseGtrees, int numGtrees, int **index, int ***listEquals, int *sizeListEquals, int **listMerged, int prevSnap, int currSnap)
{
  int gal = 0;
  int count = 0;
  int indexGal = -1;
  int indexDescGal = -1;
  int *thisIndex = NULL;
  int listContainingIndexProg = 0;
  int listContainingIndexDesc = 0;

  outgtree_t *thisGtree = NULL;
  outgal_t thisGal;
                  
  for(int tree=0; tree<numGtrees; tree++)
  {
    gal = 0;
    thisGtree = theseGtrees[tree];
    thisIndex = index[tree];    // 2D array with dimenions [numTree][numGal]
    count = get_max(thisIndex, thisGtree->numGal);
    count++;                 
        
    /* find galaxies that is at prevSnap */
    while(thisGtree->galaxies[gal].snapnumber < prevSnap + 1)
      gal++;
    
    /* go through galaxies until reaching currSnap+1 */
    while(thisGtree->galaxies[gal].snapnumber < currSnap + 1 && gal < thisGtree->numGal)
    {
      thisGal = thisGtree->galaxies[gal];
      
      /* if galaxies is a new branch */
      if(thisIndex[gal] == -1)
      {
        thisIndex[gal] = count;
        count++;
        add_index(&(listEquals[tree]), &(sizeListEquals[tree]), thisIndex[gal]);
      }
      /* check if galaxy needs to be merged */
      else
      {
        while(is_index_in_toBeMergedList(thisIndex[gal], &(listMerged[tree])) == 1)
        {
          /* merge */
          indexGal = get_index_toBeMergedList(thisIndex[gal], &(listMerged[tree]));
          indexDescGal = thisIndex[gal];
          listContainingIndexProg = search_index(indexGal, listEquals[tree]);
          listContainingIndexDesc = search_index(indexDescGal, listEquals[tree]);
          merge_lists(&(listEquals[tree]), &(sizeListEquals[tree]), listContainingIndexProg, listContainingIndexDesc);
        }
      }
      
      if(thisGal.localDescID < 0)
        break;

      /* if descendent galaxy has NOT been indexed before */
      if(thisIndex[thisGal.localDescID] == -1)
      {
        thisIndex[thisGal.localDescID] = thisIndex[gal];
      }
      /* if descendent galaxy has been indexed before */
      else
      {
        /* write to list to be merged */
        write_in_toBeMergedList(thisIndex[gal], thisIndex[thisGal.localDescID], &(listMerged[tree]));
      }
      
      gal++;
    }
  }
}