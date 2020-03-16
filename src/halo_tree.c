#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "halo_tree.h"


/*-------------------------------------------------------*/
/* FUNCTIONS ON HALO_T                                   */
/*-------------------------------------------------------*/

halo_t *initHalo()
{
  halo_t *newHalo = malloc(sizeof(halo_t));
  if(newHalo == NULL)
  {
    fprintf(stderr, "Could not allocate halo_t\n");
    exit(EXIT_FAILURE);
  }

  newHalo->ID = -1;
  newHalo->descID = -1;
  newHalo->localID = -1;
  newHalo->localDescID = -1;
  newHalo->numProg = 0;
  
  newHalo->snapnumber = -1;
  newHalo->scalefactor = -1.;
  newHalo->descScalefactor = -1.;
  
  for(int i=0; i<3; i++)
  {
    newHalo->pos[i] = 0.;
    newHalo->vel[i] = 0.;
  }
  newHalo->Mvir = 0.;
  newHalo->Rvir = 0.;
  newHalo->velDisp = 0.;
  newHalo->velMax = 0.;
  newHalo->halfmassRadius = 0.;
  newHalo->spin = 0.;
  newHalo->spinBullock = 0.;
  newHalo->scalefactorLastMajorMerger = 0.;
  
  return newHalo;
}

void deallocate_halo(halo_t *thisHalo)
{
  if(thisHalo != NULL) free(thisHalo);
}

/*-------------------------------------------------------*/
/* FUNCTIONS ON TREE_T                                   */
/*-------------------------------------------------------*/

tree_t *initTree(int Nhalos)
{
  tree_t *newTree = malloc(sizeof(tree_t));
   if(newTree == NULL)
  {
    fprintf(stderr, "Could not allocate tree_t.\n");
    exit(EXIT_FAILURE);
  }
  
  newTree->numHalos = Nhalos;
  newTree->halos = malloc(sizeof(halo_t) * Nhalos);
  if(newTree == NULL)
  {
    fprintf(stderr, "Could not allocate halos (Nhalos * halo_t) in tree.\n");
    exit(EXIT_FAILURE);
  }
  
  return newTree;
}

void reallocTree(tree_t **thisTree, int Nhalos)
{  
  halo_t *newHalos = NULL;
  tree_t *newTree = NULL;
  
  if((*thisTree)->halos != NULL)
  {
    newHalos = realloc((*thisTree)->halos, sizeof(halo_t) * Nhalos);
    (*thisTree)->numHalos = Nhalos;
    (*thisTree)->halos = newHalos;
  }
  else
  {
    newTree = initTree(Nhalos);
    *thisTree = newTree;
  }
}

void deallocate_tree(tree_t *thisTree)
{
  if(thisTree->halos != NULL) free(thisTree->halos);
  thisTree->halos = NULL;
  if(thisTree != NULL) free(thisTree);
  thisTree = NULL;
}

void deallocate_treeList(tree_t **theseTrees, int numTrees)
{
  for(int tree=0; tree<numTrees; tree++)
  {
    deallocate_tree(theseTrees[tree]);
  }
  free(theseTrees);
}
