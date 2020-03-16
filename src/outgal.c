#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "outgal.h"

/*-------------------------------------------------------*/
/* FUNCTIONS ON OUTGALSNAP_T                             */
/*-------------------------------------------------------*/

outgalsnap_t *initOutGalSnap()
{
  outgalsnap_t *newOutGal = malloc(sizeof(outgalsnap_t));
  if(newOutGal == NULL)
  {
    fprintf(stderr, "Could not allocate outgalsnap_t\n");
    exit(EXIT_FAILURE);
  }

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

  for(int i=0; i<73; i++)
    newOutGal->stellarmasshistory[i] = 0.;

  return newOutGal;
}

outgalsnap_t *initOutGalList(int numGal)
{
  outgalsnap_t *newOutGalList = malloc(sizeof(outgalsnap_t) * numGal);
  if(newOutGalList == NULL)
  {
    fprintf(stderr, "Could not allocate newOutGalList\n");
    exit(EXIT_FAILURE);
  }
  
  for(int gal=0; gal<numGal; gal++)
  {    
    newOutGalList[gal].scalefactor = -1.;

    for(int i=0; i<3; i++)
    {
      newOutGalList[gal].pos[i] = 0.;
      newOutGalList[gal].vel[i] = 0.;
    }
    newOutGalList[gal].Mvir = 0.;
    newOutGalList[gal].Mvir_prog = 0.;
    newOutGalList[gal].Rvir = 0.;
    newOutGalList[gal].velDisp = 0.;
    newOutGalList[gal].velMax = 0.;
    newOutGalList[gal].spin = 0.;
    newOutGalList[gal].scalefactorLastMajorMerger = 0.;
    newOutGalList[gal].Mgas = 0.;
    newOutGalList[gal].MgasIni = 0.;
    newOutGalList[gal].Mstar = 0.;
    newOutGalList[gal].feff = 0.;
    newOutGalList[gal].fg = 0.;
    newOutGalList[gal].zreion = 0.;
    newOutGalList[gal].photHI_bg = 0.;

    for(int i=0; i<73; i++)
      newOutGalList[gal].stellarmasshistory[i] = 0.;

  }
  
  return newOutGalList;
}

void deallocate_outgalsnap(outgalsnap_t *thisOutGal)
{
  if(thisOutGal != NULL) free(thisOutGal);
}

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
    fprintf(stderr, "Could not allocate outgtree_t.\n");
    exit(EXIT_FAILURE);
  }
  
  newOutGtree->numGal = numGal;
  newOutGtree->outgalaxies = malloc(sizeof(outgal_t) * numGal);
  if(newOutGtree == NULL)
  {
    fprintf(stderr, "Could not allocate galaxies (numGal * outgal_t) in tree.\n");
    exit(EXIT_FAILURE);
  }
  
  return newOutGtree;
}

void deallocate_outgtree(outgtree_t *thisOutGtree)
{
  if(thisOutGtree->outgalaxies != NULL) free(thisOutGtree->outgalaxies);
  thisOutGtree->outgalaxies = NULL;
  if(thisOutGtree != NULL) free(thisOutGtree);
  thisOutGtree = NULL;
}
