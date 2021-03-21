#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "gal_gtree.h"
#include "cutandresort.h"

/*
Cut all the trees at outSnap and resort them. For each trees, associate a unique index to each galaxy at the beginning of each branch and propagate it. Keep track of all the index that are contained in the same branch. When two branches merge, merge the lists containing the two branch index. At the end, walk the original tree again, sorting all the galaxies in the new trees depending on their index. 
*/
double cutandresort(gtree_t **theseGtrees, int numGtrees, int outSnap, gtree_t ***newGtrees)
{
    int numberOfRoots = 0;
    int j = 0;    
    int count = 0;
    int sizeListEquals = 0;
    int **listEquals = NULL;
    int *index = NULL;
    int treeToPutThisGalIn;
    int listContainingIndexProg;
    int listContainingIndexDesc;
    int numberNewGtrees = 0;
    gtree_t **theseNewGtrees = NULL;
    gal_t thisGal;
    gtree_t **new_trees = NULL;
    
    for(int i=0; i<numGtrees; i++)
    {        
        index = malloc(sizeof(int) * theseGtrees[i]->numGal);
        for(int h=0; h<theseGtrees[i]->numGal; h++)
            index[h] = -1;
        numberOfRoots = get_number_of_root(theseGtrees[i], outSnap);
        if(numberOfRoots != 0)
        {
            j = 0;
            count = 0;
            sizeListEquals = 0;
            
            while(theseGtrees[i]->galaxies[j].snapnumber < (outSnap+1))
            {
                thisGal = theseGtrees[i]->galaxies[j];

                if(index[j] == -1)
                {
                    index[j] = count;
                    count++;
                    add_index(&listEquals, &sizeListEquals, index[j]);
                }

                if((index[thisGal.localDescID] == -1) || (thisGal.snapnumber == outSnap))
                    index[thisGal.localDescID] = index[j];
                else
                {
                    listContainingIndexProg = search_branch(index[j], listEquals);
                    listContainingIndexDesc = search_branch(index[thisGal.localDescID], listEquals);
                    merge_lists(&listEquals, &sizeListEquals, listContainingIndexProg, listContainingIndexDesc);
                }
                
                if(thisGal.snapnumber == outSnap)
                    (theseGtrees[i]->galaxies+j)->localDescID = -1;
                
                j++;

            }
            
            if(numberOfRoots != sizeListEquals)
            {
                printf("for tree %d : different size in test.c line 54 !\n number of roots = %d and size of list of equivalent index = %d\n", i, numberOfRoots, sizeListEquals);
                exit(0);
            }
            
            new_trees = malloc(sizeof(gtree_t*) * sizeListEquals);
            for(int t=0; t<sizeListEquals; t++)
            {
                new_trees[t] = initGtree(1); // have to initialize at 1 gal
                new_trees[t]->numGal = -1;
            }

            for(int l=0; l<j; l++)
            {
                treeToPutThisGalIn = search_branch(index[l], listEquals);
                add_gal_to_tree(&new_trees, treeToPutThisGalIn, theseGtrees[i]->galaxies[l]);
            }    

            reattribute_localID(new_trees, sizeListEquals);
            add_trees(&theseNewGtrees, &numberNewGtrees, new_trees, sizeListEquals);        
            
            for(int m=0; m<sizeListEquals; m++)
            {
                free(new_trees[m]);
                free(listEquals[m]);
            }
            free(listEquals);
            listEquals = NULL;
            
            free(new_trees);
            new_trees = NULL;
        }
        free(index);
    }

    *newGtrees = theseNewGtrees;
    return numberNewGtrees;
}

int search_index(int *arrayToSearch, int value)
{
    int count = 0;
    while(arrayToSearch[count] != value)
        count++;

    return count;
}      

void reattribute_localID(gtree_t **thisGtrees, int numberTrees)
{
    int *tmpArray = NULL;

    for(int i=0; i<numberTrees; i++)
    {
        int *ptr_tmp = realloc(tmpArray, sizeof(int) * thisGtrees[i]->numGal);
        if(ptr_tmp == NULL)
        {
            fprintf(stderr, "realloc failed, exiting !\n");
            free(tmpArray);
            exit(EXIT_FAILURE);
        }
        else
            tmpArray = ptr_tmp;

        for(int j=0; j<thisGtrees[i]->numGal; j++)
            tmpArray[j] = thisGtrees[i]->galaxies[j].localID;
        
        for(int j=0; j<thisGtrees[i]->numGal-1; j++)
        {
            (thisGtrees[i]->galaxies+j)->localID = search_index(tmpArray, thisGtrees[i]->galaxies[j].localID);
            (thisGtrees[i]->galaxies+j)->localDescID = search_index(tmpArray, thisGtrees[i]->galaxies[j].localDescID);
        }

        (thisGtrees[i]->galaxies+thisGtrees[i]->numGal-1)->localID = search_index(tmpArray, thisGtrees[i]->galaxies[thisGtrees[i]->numGal-1].localID);
        (thisGtrees[i]->galaxies+thisGtrees[i]->numGal-1)->localDescID = -1;
       }
       
       free(tmpArray);
}

int get_number_of_root(gtree_t *thisGtree, int outSnap)
{
    int numberOfRoot = 0;
    int i = 0;
    while(thisGtree->galaxies[i].snapnumber < (outSnap+1))
    {
        if(thisGtree->galaxies[i].snapnumber == outSnap)
            numberOfRoot++;
        i++;
    }
            
    return numberOfRoot;
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

int search_branch(int index, int **listEquals)
{
    int wantedList = 0;

    while(is_not_in_this_list(index, listEquals[wantedList]))
        wantedList++;
    
    return wantedList;
}

int is_not_in_this_list(int index, int *thisList)
{
    int isNot = 1;
    for(int i=1; i<thisList[0]+1; i++)
        if(index == thisList[i])
            isNot = 0;
    
    return isNot;
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

void add_gal_to_tree(gtree_t ***theseGtree, int position, gal_t thisGal)
{
    gtree_t *thisGtree = (*theseGtree)[position];
    
    int i = thisGtree->numGal;
    if(i==-1)
        i++;    
    reallocGtree(&thisGtree, i + 1);

    
    thisGtree->galaxies[i] = thisGal; //careful, both stellarmasshistory point to the same data
}

void add_trees(gtree_t ***theseNewGtrees, int *numberNewGtrees, gtree_t **new_trees, int sizeListEquals)
{
    gtree_t **localNewGtree = *theseNewGtrees;

    localNewGtree = realloc(*theseNewGtrees, sizeof(gtree_t*) * (*numberNewGtrees + sizeListEquals));
    if(localNewGtree == NULL)
    {
        printf("Could not realloc NewGtree ! Exiting ! \n");
        exit(0);
    }
    for(int i=0; i<sizeListEquals; i++)
    {
        localNewGtree[i + *numberNewGtrees] = new_trees[i];
        new_trees[i] = NULL; // so we can deallocate new_trees without deallocating localNewGtree 
    }
        
    *numberNewGtrees += sizeListEquals;
    *theseNewGtrees = localNewGtree;
}

