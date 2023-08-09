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
  
  newOutGal->MgasIni = 0.;
#ifndef FIRST
  newOutGal->fracMgasMer = 0.;
#endif
#if defined WITHMETALS
  newOutGal->MgasNew = 0.;
  newOutGal->MgasEj = 0.;
#endif
  newOutGal->Mgas = 0.;
  newOutGal->Mstar = 0.;
  
#ifndef FIRST
  newOutGal->fesc = 0.;
  newOutGal->Nion = 0.;
  newOutGal->fej = 0.;
#endif
  newOutGal->feff = 0.;
  newOutGal->fg = 0.;
  newOutGal->zreion = 0.;  
  newOutGal->photHI_bg = 0.;  

  for(int i=0; i<OUTPUTSNAPNUMBER; i++)
    newOutGal->stellarmasshistory[i] = 0.;

#if defined WITHMETALS
  for(int i=0; i<3; i++)
  {
    newOutGal->Mmetal[i] = 0.;
    newOutGal->MmetalIni[i] = 0.;
    newOutGal->MmetalNew[i] = 0.;
    newOutGal->fracMmetalMer[i] = 0.;
    newOutGal->MmetalEj[i] = 0.;
  }

  for(int i=0; i<OUTPUTSNAPNUMBER; i++)
    newOutGal->metalmasshistory[i] = 0.;
  
  newOutGal->Mdust = 0.;
#endif

  return newOutGal;
}

void initOutGalSnapWithoutAllocation(outgalsnap_t *newOutGal)
{
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
  
  newOutGal->MgasIni = 0.;
#ifndef FIRST
  newOutGal->fracMgasMer = 0.;
#endif
#if defined WITHMETALS
  newOutGal->MgasNew = 0.;
  newOutGal->MgasEj = 0.;
#endif
  newOutGal->Mgas = 0.;
  newOutGal->Mstar = 0.;
  
#ifndef FIRST
  newOutGal->fesc = 0.;
  newOutGal->Nion = 0.;
  newOutGal->fej = 0.;
#endif
  newOutGal->feff = 0.;
  newOutGal->fg = 0.;
  newOutGal->zreion = 0.;  
  newOutGal->photHI_bg = 0.;  

  for(int i=0; i<OUTPUTSNAPNUMBER; i++)
    newOutGal->stellarmasshistory[i] = 0.;

#if defined WITHMETALS
  for(int i=0; i<3; i++)
  {
    newOutGal->Mmetal[i] = 0.;
    newOutGal->MmetalIni[i] = 0.;
    newOutGal->MmetalNew[i] = 0.;
    newOutGal->fracMmetalMer[i] = 0.;
    newOutGal->MmetalEj[i] = 0.;
  }

  for(int i=0; i<OUTPUTSNAPNUMBER; i++)
    newOutGal->metalmasshistory[i] = 0.;
  
  newOutGal->Mdust = 0.;
#endif
}

outgalsnap_t *initOutGalList(int numGal)
{
  outgalsnap_t *newOutGalList = malloc(sizeof(outgalsnap_t) * numGal);
  if(newOutGalList == NULL)
  {
    fprintf(stderr, "Could not allocate newOutGalList\n");
    exit(EXIT_FAILURE);
  }
  memset(newOutGalList, 0, sizeof(outgalsnap_t) * numGal);
  
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
    
    newOutGalList[gal].MgasIni = 0.;
#ifndef FIRST
    newOutGalList[gal].fracMgasMer = 0.;
#endif
#if defined WITHMETALS
    newOutGalList[gal].MgasNew = 0.;
    newOutGalList[gal].MgasEj = 0.;
#endif
    newOutGalList[gal].Mgas = 0.;
    newOutGalList[gal].Mstar = 0.;
    
#ifndef FIRST
    newOutGalList[gal].fesc = 0.;
    newOutGalList[gal].Nion = 0.;
    newOutGalList[gal].fej = 0.;
#endif
    newOutGalList[gal].feff = 0.;
    newOutGalList[gal].fg = 0.;
    newOutGalList[gal].zreion = 0.;
    newOutGalList[gal].photHI_bg = 0.;

    for(int i=0; i<OUTPUTSNAPNUMBER; i++)
      newOutGalList[gal].stellarmasshistory[i] = 0.;

#ifdef WITHMETALS
    for(int i=0; i<3; i++)
    {
      newOutGalList[gal].Mmetal[i] = 0.;
      newOutGalList[gal].MmetalIni[i] = 0.;
      newOutGalList[gal].MmetalNew[i] = 0.;
      newOutGalList[gal].fracMmetalMer[i] = 0.;
      newOutGalList[gal].MmetalEj[i] = 0.;
    }

    for(int i=0; i<OUTPUTSNAPNUMBER; i++)
      newOutGalList[gal].metalmasshistory[i] = 0.;
    
    newOutGalList[gal].Mdust = 0.;
#endif
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
  
#if defined WITHROCKSTARID
  newOutGal->ID = 0;
  newOutGal->descID = 0;
#endif
  
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
  
  newOutGal->MgasIni = 0.;
#ifndef FIRST
  newOutGal->fracMgasMer = 0.;
#endif
#if defined WITHMETALS
  newOutGal->MgasNew = 0.;
  newOutGal->MgasEj = 0.;
#endif
  newOutGal->Mgas = 0.;
  newOutGal->Mstar = 0.;
  
#ifndef FIRST
  newOutGal->fesc = 0.;
  newOutGal->Nion = 0.;
  newOutGal->fej = 0.;
#endif
  newOutGal->feff = 0.;
  newOutGal->fg = 0.;
  newOutGal->zreion = 0.;
  newOutGal->photHI_bg = 0.;

#if defined WITHMETALS
  for(int i=0; i<3; i++)
  {
    newOutGal->Mmetal[i] = 0.;
    newOutGal->MmetalIni[i] = 0.;
    newOutGal->MmetalNew[i] = 0.;
    newOutGal->fracMmetalMer[i] = 0.;
    newOutGal->MmetalEj[i] = 0.;
  }
  
  newOutGal->Mdust = 0.;
#endif
  
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
  if(newOutGtree->outgalaxies == NULL)
  {
    fprintf(stderr, "Could not allocate galaxies (numGal * outgal_t) in tree.\n");
    exit(EXIT_FAILURE);
  }
  memset(newOutGtree->outgalaxies, 0, sizeof(outgal_t) * numGal);
  
  return newOutGtree;
}

void deallocate_outgtree(outgtree_t *thisOutGtree)
{
  if(thisOutGtree->outgalaxies != NULL) free(thisOutGtree->outgalaxies);
  thisOutGtree->outgalaxies = NULL;
  if(thisOutGtree != NULL) free(thisOutGtree);
  thisOutGtree = NULL;
}
