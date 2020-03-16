#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "gal_gtree.h"


/*-------------------------------------------------------*/
/* FUNCTIONS ON GAL_T                                    */
/*-------------------------------------------------------*/

gal_t *initGal()
{
  gal_t *newGal = malloc(sizeof(gal_t));
  if(newGal == NULL)
  {
    fprintf(stderr, "Could not allocate gal_t\n");
    exit(EXIT_FAILURE);
  }

  newGal->localID = -1;
  newGal->localDescID = -1;
  newGal->numProg = 0;
  
  newGal->snapnumber = -1;
  newGal->scalefactor = -1.;
  newGal->descScalefactor = -1.;
  
  for(int i=0; i<3; i++)
  {
    newGal->pos[i] = 0.;
    newGal->vel[i] = 0.;
  }
  newGal->Mvir = 0.;
  newGal->Mvir_prog = 0.;
  newGal->Rvir = 0.;
  newGal->velDisp = 0.;
  newGal->velMax = 0.;
  newGal->spin = 0.;
  newGal->scalefactorLastMajorMerger = 0.;
  newGal->Mgas = 0.;
  newGal->Mstar = 0.;
  newGal->MgasIni = 0.;
  
  newGal->feff = 0.;
  newGal->fg = 0.;
  newGal->photHI_bg = 0.;
  newGal->zreion = 0.;
  
  newGal->stellarmasshistory = NULL;

  printf("pointer is %p\n", newGal->stellarmasshistory);
  return newGal;
}

void deallocate_gal(gal_t *thisGal)
{
  if(thisGal->stellarmasshistory != NULL) free(thisGal->stellarmasshistory);
  if(thisGal != NULL) free(thisGal);
}

/*-------------------------------------------------------*/
/* FUNCTIONS ON GTREE_T                                  */
/*-------------------------------------------------------*/

gtree_t *initGtree(int Ngal)
{
  gtree_t *newGtree = malloc(sizeof(gtree_t));
  if(newGtree == NULL)
  {
    fprintf(stderr, "Could not allocate gtree_t.\n");
    exit(EXIT_FAILURE);
  }
  
  newGtree->numGal = Ngal;
  newGtree->walker = 0;
  newGtree->galaxies = malloc(sizeof(gal_t) * Ngal);
  if(newGtree == NULL)
  {
    fprintf(stderr, "Could not allocate galaxies (Ngal * gal_t) in tree.\n");
    exit(EXIT_FAILURE);
  }
  
  return newGtree;
}

void reallocGtree(gtree_t **thisGtree, int Ngal)
{  
  gal_t *newGals = NULL;
  gtree_t *newGtree = NULL;
  
  if((*thisGtree)->galaxies != NULL)
  {
    newGals = realloc((*thisGtree)->galaxies, sizeof(gal_t) * Ngal);
    (*thisGtree)->numGal = Ngal;
    (*thisGtree)->galaxies = newGals;
  }
  else
  {
    newGtree = initGtree(Ngal);
    *thisGtree = newGtree;
  }
}

void deallocate_gtree(gtree_t *thisGtree)
{
  if(thisGtree->galaxies != NULL) free(thisGtree->galaxies);
  thisGtree->galaxies = NULL;
  if(thisGtree != NULL) free(thisGtree);
  thisGtree = NULL;
}

void deallocate_gtreeList(gtree_t **theseGtrees, int numGtrees)
{
  for(int tree=0; tree<numGtrees; tree++)
  {
    deallocate_gtree(theseGtrees[tree]);
  }
  free(theseGtrees);
}
