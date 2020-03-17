#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "outgal.h"
#include "tree_correct.h"

void correct_localDescID_in_trees(outgtree_t **theseTrees, int32_t numTrees, int endSnap)
{
  outgtree_t *thisTree = NULL;
  outgal_t *thisGal = NULL;
  int32_t *oldLocalIDs = NULL;
  
  for(int tree=0; tree<numTrees; tree++)
  {
    thisTree = theseTrees[tree];
    
    if(thisTree->galaxies[thisTree->numGal-1].localID >= thisTree->numGal)
    {  
      oldLocalIDs = allocate_array_int32_t(thisTree->numGal, "oldLocalIDs");
      
      for(int gal=0; gal<thisTree->numGal; gal++)
      {
        thisGal = &(thisTree->galaxies[gal]);
        
        if(thisGal->localID != gal)
        {
          oldLocalIDs[gal] = thisGal->localID;
          thisGal->localID = gal;
        }
      }
      
      for(int gal=0; gal<thisTree->numGal; gal++)
      {
        thisGal = &(thisTree->galaxies[gal]);
        
        /* find descID in oldLocalIDs */
        for(int i=0; i<thisTree->numGal; i++)
        {
          if(thisGal->localDescID == oldLocalIDs[i])
          {
            thisGal->localDescID = i;
          }
        }
      }
      free(oldLocalIDs);
    }
      
    for(int gal=0; gal<thisTree->numGal; gal++)
    {
      thisGal = &(thisTree->galaxies[gal]);
      if(thisGal->snapnumber == endSnap && thisGal->localDescID >=0)
      {
        thisGal->localDescID = -1;
      }
    }
  }
}