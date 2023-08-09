#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "halo_tree.h"
#include "gal_gtree.h"
#include "copy_halo_gal.h"

void copy_halo_to_gal(halo_t *thisHalo, gal_t *thisGal)
{
#if defined WITHROCKSTARID
  thisGal->ID = thisHalo->ID;
  thisGal->descID = thisHalo->descID;
#endif
  
  thisGal->localID = thisHalo->localID;
  thisGal->localDescID = thisHalo->localDescID;
  thisGal->numProg = thisHalo->numProg;
  
  thisGal->snapnumber = thisHalo->snapnumber;
  thisGal->scalefactor = thisHalo->scalefactor;
  thisGal->descScalefactor = thisHalo->descScalefactor;
  
  for(int i=0; i<3; i++)
  {
    thisGal->pos[i] = thisHalo->pos[i];
    thisGal->vel[i] = thisHalo->vel[i];
  }
  thisGal->Mvir = thisHalo->Mvir;
  thisGal->Mvir_prog = 0.;
  thisGal->Rvir = thisHalo->Rvir;
  thisGal->velDisp = thisHalo->velDisp;
  thisGal->velMax = thisHalo->velMax;
  thisGal->spin = thisHalo->spin;
  thisGal->scalefactorLastMajorMerger = thisHalo->scalefactorLastMajorMerger;
  
  thisGal->MgasIni = 0.;
  thisGal->fracMgasMer = 0.;
#if defined WITHMETALS
  thisGal->MgasNew = 0.;
  thisGal->MgasEj = 0.;
#endif
  thisGal->Mgas = 0.;
  thisGal->Mstar = 0.;
  
  thisGal->fesc = 0.;
  thisGal->Nion = 0;
  thisGal->fej = 0.;
  thisGal->feff = 0.;
  thisGal->fg = 0.;
  thisGal->photHI_bg = 0.;
  thisGal->zreion = 0.;
  
  thisGal->stellarmasshistory = NULL;

#if defined WITHMETALS
  for(int i=0; i<3; i++)
  {
    thisGal->Mmetal[i] = 0.; 
    thisGal->MmetalIni[i] = 0.;
    thisGal->fracMmetalMer[i] = 0.;
    thisGal->MmetalNew[i] = 0.;
    thisGal->MmetalEj[i] = 0.;
    thisGal->igmMetallicity[i] = 0.;
  }
  thisGal->metalmasshistory = NULL;
  
  thisGal->Mdust = 0.;
  thisGal->MdustEj = 0.;
  thisGal->igmDustFraction = 0.;
#endif
}

gal_t *gal_from_halo(halo_t *thisHalo)
{
  gal_t *newGal = initGal();
  
  copy_halo_to_gal(thisHalo, newGal);
  
  return newGal;
}

gtree_t *gtree_from_tree(tree_t *thisTree)
{
  halo_t *thisHalo = NULL;
  gal_t *thisGal = NULL;
  int numGal = thisTree->numHalos;
  gtree_t *newGtree = initGtree(numGal);
 
  for(int gal=0; gal<numGal; gal++)
  {
    thisHalo = &(thisTree->halos[gal]);
    thisGal = &(newGtree->galaxies[gal]);
    copy_halo_to_gal(thisHalo, thisGal);
  }
  
  return newGtree;
}

int32_t gtrees_from_trees(tree_t ***thisTreeList, int numTrees, gtree_t ***thisGtreeList)
{
  tree_t *thisTree = NULL;
  tree_t **theseTrees = *thisTreeList;
  gtree_t **theseGtrees = *thisGtreeList;
  int32_t numGtrees = numTrees;
    
  theseGtrees = malloc(sizeof(gtree_t) * numGtrees);
  for(int gtree=0; gtree<numGtrees; gtree++)
  {
    thisTree = theseTrees[gtree];
    theseGtrees[gtree] = gtree_from_tree(thisTree);
    deallocate_tree(thisTree);
  }

  *thisGtreeList = theseGtrees;
    
  return numGtrees;
}
