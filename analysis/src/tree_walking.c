#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "utils.h"
#include "dconfObj.h"
#include "outgal.h"

#include "tree_walking.h"

int getmax(int *listOfInt, int length)
{
        int max = -10;
        for(int i=0; i<length; i++)
                if(listOfInt[i] > max)
                        max = listOfInt[i];
        
        return max;
}

void printlistequals(int **listEquals, int numberlist)
{
        for(int i=0; i<numberlist; i++)
        {
                printf("[");
                for(int j=1; j<listEquals[i][0]+1; j++)
                        printf(" %d ",listEquals[i][j]);
                printf("]");
        }
        printf("\n");
}

void getListEquals(outgtree_t **theseGtrees, int numGtrees, int **index, int ***listEquals, int *sizeListEquals, int prevSnap, int currentSnap)
{
        int galCount;
        int count;
        int *thisIndex;
        int list_containing_index_prog;
        int list_containing_index_desc;
        
        
        outgtree_t *thisGtree = NULL;
        outgal_t thisGal;
        
        
        if(prevSnap != 0)
                        prevSnap++; //want to include the very first snap but not in the subsequent loops
                        
        for(int i=0; i<numGtrees; i++)
        {

                galCount = 0;
                thisGtree = theseGtrees[i];
                thisIndex = index[i];
                count = getmax(thisIndex, thisGtree->numGal);
                count++;                 
                
                while(thisGtree->galaxies[galCount].snapnumber < prevSnap)
                        galCount++;

                if(i==0)
                        printf("galcountstart = %d\n", galCount);
                while(thisGtree->galaxies[galCount].snapnumber < currentSnap + 1)
                {
                        thisGal = thisGtree->galaxies[galCount];
                        if(thisIndex[galCount] == -1)
                        {
                                thisIndex[galCount] = count;
                                count++;
                                add_index(&(listEquals[i]), &(sizeListEquals[i]), thisIndex[galCount]);
                        }
                        if(thisIndex[thisGal.localDescID] == -1)
                                thisIndex[thisGal.localDescID] = thisIndex[galCount];
                        else
                        {                                       
                                list_containing_index_prog = search_index(thisIndex[galCount], listEquals[i]);
                                list_containing_index_desc = search_index(thisIndex[thisGal.localDescID], listEquals[i]);
                                mergeLists(&(listEquals[i]), &(sizeListEquals[i]), list_containing_index_prog, list_containing_index_desc);
                        }

                        galCount++;
                }
                if(i==0)
                        printf("galcountfinished = %d\n", galCount);
        }
        //exit(0);
}

void add_index(int ***listEquals, int *sizeListEquals, int index)
{
        (*sizeListEquals)++;
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

        while(isNotInThisList(index, listEquals[wantedList]))
                wantedList++;
        
        return wantedList;
}

int isNotInThisList(int index, int *thisList)
{
        int isNot = 1;
        for(int i=1; i<thisList[0]+1; i++)
                if(index == thisList[i])
                        isNot = 0;
        
        return isNot;
}

// merge two lists in the listEquals, free the second list and move the list with higher index accordingly
void mergeLists(int ***listEquals, int *sizeListEquals, int listOne, int listTwo)
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
