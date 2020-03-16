#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"

#include "dconfObj.h"
#include "halo_tree.h"
#include "gal_gtree.h"
#include "copy_halo_gal.h"
#include "read_trees.h"
#include "outgal.h"
#include "output.h"
#include "redshift_list.h"

#include "domain.h"
#include "nion.h"
#include "nion_gal.h"
#include "ion_emissivity.h"
#include "run.h"

#include "cifog/confObj.h"
#include "cifog/grid.h"
#include "cifog/photion_background.h"
#include "cifog/recombination.h"
#include "cifog/cifog_delphi.h"

#include "check_cifogParam.h"
#include "photion_gal.h"
#include "evol_gal.h"
#include "starformation_SNfeedback.h"

#include "cutandresort.h"
#define MAXGAL 100000
#define MAXGAL_INC 10000

void run_delphi(dconfObj_t simParam, int thisRank, int size)
{
  /*---------------------------------------------------------------------*/
  /* INITIALIZATION & STARTUP */
  /*---------------------------------------------------------------------*/

  /* INITILAIZE DOMAIN */
  domain_t *domain = initDomain(simParam->gridsize, thisRank, size);

  /* LOADING TREES */
  tree_t **theseTrees = NULL;
  int32_t numTrees = read_trees_in_all_files(simParam->fileName, simParam->numFiles, &theseTrees, thisRank, size);

  /* COPYING TREE TO GALAXY STRUCT & DELETE HALO TREES */
  gtree_t **theseGtrees = NULL;
  int32_t numGtrees = gtrees_from_trees(&theseTrees, numTrees, &theseGtrees);
  free(theseTrees);

  /* SIMULATION */
  redshiftlist_t *thisRedshiftList = read_redshift_snap_list(simParam->redshiftFile);
  int startSnap = simParam->startSnap;
  int endSnap = simParam->endSnap;
  int deltaSnap = simParam->deltaSnap;
  int minSnap = 0;
  int maxSnap = 0;
  simParam->timeSnaps = get_times_from_redshifts(thisRedshiftList, endSnap, simParam->hubble_h, simParam->omega_m, simParam->omega_l);
  
  if(simParam->delayedSNfeedback == 1)
    simParam->SNenergy = get_SNenergy_delayed(endSnap, simParam->timeSnaps);
  
  if(simParam->reion == 1 && simParam->sps_model > 10)
    simParam->corrFactor_nion = get_corrFactor_nion(simParam);
  
  /* WRITING OUTPUT */
  int *numGalSnap = allocate_array_int(endSnap, "numGalSnap");
  outgalsnap_t *outGalList = NULL;
  int outSnapCounter = 0;

  /* REIONIZATION PARAMETERS */
  int memory_intensive = 1;
  nion_t *thisNionList = NULL;
  confObj_t cifogParam = NULL;
  grid_t *grid = NULL;
  fftw_complex *cifogNion = NULL;
  fftw_complex *cifogNionHeI = NULL;
  fftw_complex *cifogNionHeII = NULL;
  fftw_complex *photHI = NULL;
  fftw_complex *XHII = NULL;
  fftw_complex *zreion = NULL;
  integral_table_t *integralTable = NULL;
  photIonlist_t *photIonBgList = NULL;
  if(simParam->reion == 1 && simParam->reion_model == 2)
  {
    cifog_init(simParam->cifogIniFile, &cifogParam, &grid,  &integralTable, &photIonBgList, thisRank);
    check_cifogParam(simParam, cifogParam, thisRank);
  }
  
  /* SET WALKER TO CORRECT STARTING POSITION IN TREE */
  set_walker_to_startSnap(numGtrees, &theseGtrees, startSnap);
  
  int numHalos = 0;     // currently only for printing output needed
  
  /*---------------------------------------------------------------------*/
  /* RUN DELPHI ACROSS SNAPSHOTS */
  /*---------------------------------------------------------------------*/
  
  /* LOOPING OVER SNAPSHOTS */
  for(int snap=startSnap; snap<endSnap; snap=snap+deltaSnap)
  {    
    if(snap == simParam->snapList[outSnapCounter])
      outGalList = initOutGalList(MAXGAL);
        
    if(simParam->reion == 1)
      thisNionList = initNion(MAXGAL);

    minSnap = snap;
    maxSnap = snap + deltaSnap;
    
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    evolve(numGtrees, &theseGtrees, minSnap, maxSnap, simParam, domain, thisNionList, photHI, XHII, zreion, simParam->snapList[outSnapCounter], &outGalList, numGalSnap);
    
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
    if(snap == simParam->snapList[outSnapCounter])
    {
      /* write galaxies in snap to file */
      write_galaxies_of_snap_to_file(simParam, snap, numGalSnap[snap], outGalList);
      free(outGalList);
      
      if(outSnapCounter < simParam->numSnapsToWrite - 1)
        outSnapCounter++;
    }

#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
        
    /* do reionization */
    if(simParam->reion == 1)
    {
      reallocNion(&thisNionList, thisNionList->numGalWritten);

      numHalos += thisNionList->numGalWritten;
      
      // print mean number of ionizing photons
      double meanNion = get_mean_nion(thisNionList);
      if(thisRank == 0) printf("\n********\nDELPHI: sum(nion) = %e\n********\n", 1.e50*meanNion);
      
      // mapping ionizing photons to grid
      cifogNion = map_galnion_to_grid(thisNionList, domain, memory_intensive);
      
      if(simParam->reion_model == 2)
      {
        // do cifog_step which should return photoionization rate grid
        double deltaRedshift = calc_deltaRedshift(thisRedshiftList, snap);
        cifog(cifogParam, snap, thisRedshiftList->redshifts[snap], deltaRedshift, grid, cifogNion, cifogNionHeI, cifogNionHeII, integralTable, photIonBgList, thisRank);
        
       // update current photoionization & XHII field
//         update_photHI_field(&photHI, grid);
//         update_XHII_field(&XHII, grid);
        update_photHI_zreion_field(&photHI, grid, thisRedshiftList->redshifts[snap], simParam->XHII_threshold);
        update_zreion_field(&zreion, grid, thisRedshiftList->redshifts[snap], simParam->XHII_threshold);
        
        char zreionFile[MAXLENGTH];
        if(thisRank == 0)
        {
          for(int i=0; i<MAXLENGTH; i++) zreionFile[i] = '\0';
          sprintf(zreionFile, "%s/zreion_%02d", simParam->outFileName, snap);
          write_grid_to_file_double_test(zreion, grid->nbins, grid->nbins, zreionFile);
        }
        char photHIFile[MAXLENGTH];
        if(thisRank == 0)
        {
          for(int i=0; i<MAXLENGTH; i++) photHIFile[i] = '\0';
          sprintf(photHIFile, "%s/photHI_zreion_%02d", simParam->outFileName, snap);
          write_grid_to_file_double_test(photHI, grid->nbins, grid->nbins, photHIFile);
        }
      }
      if(cifogNion != NULL) fftw_free(cifogNion);
      deallocate_nion(thisNionList);
      
    }
    
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
        
    /* printing */
    if(simParam->horizontalOutput==1)
    {
      int printNumGalSnap = numGalSnap[snap];
#ifdef MPI
      MPI_Allreduce(&(numGalSnap[snap]), &printNumGalSnap, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
      if(thisRank == 0) printf("snap %d: numGalSnap = %d\n", snap, printNumGalSnap);
    }
  }

  /* Write number of galaxies contained in each snap into a file */
  if(simParam->horizontalOutput==1)
  {
    write_num_galaxies_to_file(simParam, numGalSnap, endSnap, thisRank);
  }
  
  /* Write trees into #PROC file */
  if(simParam->verticalOutput == 1)
  {
    if(simParam->endSnap!=74)
    {
      gtree_t **theseNewGtrees = NULL;
      int newNumGtrees;
      
      newNumGtrees = cutandresort(theseGtrees, numGtrees, endSnap-1, &theseNewGtrees);
      write_treelist(simParam, thisRank, newNumGtrees, theseNewGtrees);
      deallocate_gtreeList(theseNewGtrees, newNumGtrees);
      printf("Rank %d : total trees = %d\n", thisRank, newNumGtrees);
    }
    else
    {
      write_treelist(simParam, thisRank, numGtrees, theseGtrees);
      printf("Rank %d : total trees = %d\n", thisRank, numGtrees);
    }
  }
  
#ifdef MPI
  int tmpNumHalos = numHalos;
  MPI_Allreduce(&tmpNumHalos, &numHalos, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if(thisRank == 0) printf("Total number of halos/galaxies = %d\n", numHalos);
    
  /*---------------------------------------------------------------------*/
  /* DEALLOCATION */
  /*---------------------------------------------------------------------*/

  // deallocate initialisation of cifog  end shut it down
  if(simParam->reion == 1)
  {
    if(simParam->reion_model == 2)
    {
      cifog_deallocate(cifogParam, grid, integralTable, photIonBgList, thisRank);
      if(photHI != NULL) fftw_free(photHI);
      if(XHII != NULL) fftw_free(XHII);
      if(zreion != NULL) fftw_free(zreion);
    }
  }

  deallocate_redshiftlist(thisRedshiftList);
  free(numGalSnap);
  deallocate_gtreeList(theseGtrees, numGtrees);
  deallocate_domain(domain);
}


void evolve(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, dconfObj_t simParam, domain_t *thisDomain, nion_t *thisNionList, fftw_complex *photHI, fftw_complex *XHII, fftw_complex *zreion, int outSnap, outgalsnap_t **outGalList, int *numGalInSnapList)
{
  gtree_t **theseGtrees = *thisGtreeList;
  gtree_t *thisGtree = NULL;
  gal_t *thisGal;
  int32_t gal = 0;
  int32_t snapGal = 0;
  int32_t numGal = 0;
  
  double inv_boxsize = simParam->inv_boxsize;
  
  int32_t numGalInSnap = 0;
  int32_t numGalInOutSnap = 0;
  int32_t count = 0;
  
  /* LOOP OVER TREES */
  for(int gtree=0; gtree<numGtrees; gtree++)
  {
    thisGtree = theseGtrees[gtree];
    numGal = thisGtree->numGal;
    gal = thisGtree->walker;
    
    if(gal < numGal) 
    {
      snapGal = thisGtree->galaxies[gal].snapnumber;

      /* LOOP OVER SNAPSHOTS */
      for(int snap=minSnap; snap<maxSnap; snap++)
      {
        if(snapGal == snap)
        {
          /* LOOP OVER GALAXIES */
          thisGal = &(thisGtree->galaxies[gal]);
          while(thisGal->snapnumber == snap)
          {
            /* get photoionization rate */
            if(simParam->reion == 1 && simParam->radfeedback == 1)
            {
              thisGal->photHI_bg = get_photHI_gal(thisGal, thisDomain->nbins, inv_boxsize, photHI);
              if(thisGal->zreion <= 0.) 
              {
//                 thisGal->zreion = get_zreion_gal(thisGal, thisDomain->nbins, inv_boxsize, XHII, simParam->XHII_threshold);
                thisGal->zreion = get_zreionField_gal(thisGal, thisDomain->nbins, inv_boxsize, zreion);
              }
            }
            
            /* evolve galaxy */
            evolve_gal(thisGal, thisGtree, simParam);
            
            if(thisGal->zreion > 0.) assert(thisGal->zreion >= 1./thisGal->scalefactor - 1.);

            /* store ionizing properties for reionization, TODO: make sure only galaxies at specified snapshot */
            if(simParam->reion == 1)
              copy_gal_to_nion(thisGal, simParam, thisDomain, 1.1, thisNionList);
            
            /* clean galaxy (free not needed memory) */
            clean_gal(thisGal, outSnap);
            
            gal++;
            numGalInSnap++;
            
            /* adding galaxies to list to write into file */
            if((*outGalList) != NULL && snap == outSnap)
            {
              copy_gal_to_outGalList(thisGal, numGalInOutSnap, *outGalList);
              numGalInOutSnap++;

              if(numGalInOutSnap == (MAXGAL + count * MAXGAL_INC) - 1)
              {
                  count++;
                  (*outGalList) = realloc((*outGalList), sizeof(outgalsnap_t)*(MAXGAL + count * MAXGAL_INC));
              }
            }
            
            /* stop when reaching the end of the tree */
            if(gal == numGal)
              break;
            
            /* go to the next galaxy */
            thisGal = &(thisGtree->galaxies[gal]);
          }
          if(gal < numGal)
            snapGal = thisGtree->galaxies[gal].snapnumber;
        }
        if(numGalInOutSnap > 0) assert(numGalInOutSnap == numGalInSnap);
        numGalInSnapList[snap] = numGalInSnap;
      }
      /* set walker, so that in next snap loop we know where to continue */
      thisGtree->walker = gal;
    }
  }
}

void set_walker_to_startSnap(int numGtrees, gtree_t ***thisGtreeList, int startSnap)
{
  gtree_t **theseGtrees = *thisGtreeList;
  gtree_t *thisGtree = NULL;
  int32_t gal = 0;
  int32_t snapGal = 0;
  int32_t numGal = 0;
      
  /* LOOP OVER TREES */
  for(int gtree=0; gtree<numGtrees; gtree++)
  {
    thisGtree = theseGtrees[gtree];
    numGal = thisGtree->numGal;
    gal = thisGtree->walker;
    snapGal = thisGtree->galaxies[gal].snapnumber;

    while(gal < numGal-1 && snapGal < startSnap)
    {
      gal++;
      snapGal = thisGtree->galaxies[gal].snapnumber;
    }
    
    thisGtree->walker = gal;
  }
}
