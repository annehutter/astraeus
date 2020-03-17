#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "outgal.h"

/*-------------------------------------------------------*/
/* FUNCTIONS ON OUTGAL_T                                 */
/*-------------------------------------------------------*/

outgal_t *initOutGalTree()
{
  outgal_t *newOutGal = malloc(sizeof(outgal_t));
  if(newOutGal == NULL)
  {
    fprintf(stderr, "Could not allocate outgalsnap_t\n");
    exit(EXIT_FAILURE);
  }
  newOutGal->localID = 0;
  newOutGal->localDescID = 0;
  newOutGal->numProg = 0;
  newOutGal->snapnumber = 0;
  newOutGal->scalefactor = -1.;
  
  for(int i=0; i<3; i++)
  {
    newOutGal->pos[i] = 0.;
    newOutGal->vel[i] = 0.;
  }
  newOutGal->Mvir = 0.;
  newOutGal->Mvir_prog = 0.;
  newOutGal->Rvir = 0.;
  newOutGal->velDisp = 0.;
  newOutGal->velMax = 0.;
  newOutGal->spin = 0.;
  newOutGal->scalefactorLastMajorMerger = 0.;
  newOutGal->Mgas = 0.;
  newOutGal->MgasIni = 0.;
  newOutGal->Mstar = 0.;
  newOutGal->feff = 0.;  
  newOutGal->fg = 0.;
  newOutGal->zreion = 0.;
  newOutGal->photHI_bg = 0.;
  
  return newOutGal;
}

void deallocate_outgaltree(outgal_t *thisOutGal)
{
  if(thisOutGal != NULL) free(thisOutGal);
}

/*-------------------------------------------------------*/
/* FUNCTIONS ON OUTGTREE_T                               */
/*-------------------------------------------------------*/

outgtree_t *initOutGtree(int numGal)
{
  outgtree_t *newOutGtree = malloc(sizeof(outgtree_t));
  if(newOutGtree == NULL)
  {
    fprintf(stderr, "Could not allocate outgtree.\n");
    exit(EXIT_FAILURE);
  }
  
  newOutGtree->numGal = numGal;
  newOutGtree->galaxies = malloc(sizeof(outgal_t) * numGal);
  if(newOutGtree == NULL)
  {
    fprintf(stderr, "Could not allocate galaxies (numGal * outgal_t) in tree.\n");
    exit(EXIT_FAILURE);
  }
  
  return newOutGtree;
}

void reallocOutGtree(outgtree_t **thisOutGtree, int numGal)
{  
  outgal_t *newGals = NULL;
  outgtree_t *newOutGtree = NULL;
  
  if((*thisOutGtree)->galaxies != NULL)
  {
    newGals = realloc((*thisOutGtree)->galaxies, sizeof(outgal_t) * numGal);
    (*thisOutGtree)->numGal = numGal;
    (*thisOutGtree)->galaxies = newGals;
  }
  else
  {
    newOutGtree = initOutGtree(numGal);
    *thisOutGtree = newOutGtree;
  }
}

void deallocate_outgtree(outgtree_t *thisOutGtree)
{
  if(thisOutGtree->galaxies != NULL) free(thisOutGtree->galaxies);
  thisOutGtree->galaxies = NULL;
  if(thisOutGtree != NULL) free(thisOutGtree);
  thisOutGtree = NULL;
}

void deallocate_outgtreeList(outgtree_t **theseOutGtrees, int numOutGtrees)
{
  for(int tree=0; tree<numOutGtrees; tree++)
  {
    deallocate_outgtree(theseOutGtrees[tree]);
  }
  free(theseOutGtrees);
}
