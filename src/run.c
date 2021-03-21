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

#include "cifog/confObj.h"
#include "cifog/grid.h"
#include "cifog/photion_background.h"
#include "cifog/recombination.h"
#include "cifog/cifog_delphi.h"

#include "check_cifogParam.h"
#include "photion_gal.h"
#include "evol_gal.h"
#include "starformation_SNfeedback.h"

#include "comm_gal_grid_struct.h"
#include "comm_gal_grid.h"

#include "cutandresort.h"

#include "run.h"

#define MAXGAL 100000
#define MAXGAL_INC 10000

/*---------------------------------------------------------------------*/
/* Core routine of astraeus: reads in data, initialises look-up tables,*/
/* loops over snapshots  to evolve galaxies and progress reionization, */
/* outputs to files, deallocation of all data */
/*---------------------------------------------------------------------*/
void run_astraeus(dconfObj_t simParam, int thisRank, int size)
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
  
  /* create lookup table for delayed SN feedback scheme */
  if(simParam->delayedSNfeedback == 1)
    simParam->SNenergy = get_SNenergy_delayed(endSnap, simParam->timeSnaps);
  
  /* create lookup table for correction factor when ionizing emissivity is distributed over time step*/
  if(simParam->reion == 1 && simParam->sps_model > 10)
    simParam->corrFactor_nion = get_corrFactor_nion(simParam);
  
  /* WRITING OUTPUT */
  int *numGalSnap = allocate_array_int(endSnap, "numGalSnap");
  outgalsnap_t *outGalList = NULL;
  int outSnapCounter = 0;

  /* REIONIZATION PARAMETERS */
  int memory_intensive = simParam->memoryIntensive;
  char *zreionName = "ZREION";
  char *photHIName = "PHOTHI";
  nion_t *thisNionList = NULL;
  confObj_t cifogParam = NULL;
  grid_t *grid = NULL;
  fftw_complex *cifogNion = NULL;
  fftw_complex *cifogNionHeI = NULL;
  fftw_complex *cifogNionHeII = NULL;
  fftw_complex *photHI = NULL;
  fftw_complex *zreion = NULL;
  integral_table_t *integralTable = NULL;
  photIonlist_t *photIonBgList = NULL;
  if(simParam->reion == 1 && simParam->reion_model == 1)
  {
    /* initialising CIFOG */
    cifog_init(simParam->cifogIniFile, &cifogParam, &grid,  &integralTable, &photIonBgList, thisRank);
    check_cifogParam(simParam, cifogParam, thisRank);
  }
  
  /* SET WALKER TO CORRECT STARTING POSITION IN TREE */
  set_walker_to_startSnap(numGtrees, &theseGtrees, startSnap);
  
  long int numHalos = 0;     // currently only for printing output needed
  
  /*---------------------------------------------------------------------*/
  /* RUN ASTRAEUS ACROSS SNAPSHOTS */
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
    
    /* evolve galaxies by <deltaSnap> snapshots */
    evolve(numGtrees, &theseGtrees, minSnap, maxSnap, simParam, domain, thisNionList, photHI, zreion, simParam->snapList[outSnapCounter], &outGalList, numGalSnap);
        
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
    if(simParam->reion == 1 && simParam->reion_model == 1)
    {
      reallocNion(&thisNionList, thisNionList->numGalWritten);

      numHalos += thisNionList->numGalWritten;
      
      /* print mean number of ionizing photons */
      double meanNion = get_mean_nion(thisNionList);
      if(thisRank == 0) printf("\n********\nDELPHI: sum(nion) = %e\n********\n", 1.e50*meanNion);
      
      /* mapping ionizing photons to grid */
      cifogNion = map_galnion_to_grid(thisNionList, domain, memory_intensive);

      /* run CIFOG for deltaRedshift which should return photoionization rate grid */
      double deltaRedshift = calc_deltaRedshift(thisRedshiftList, snap);
      cifog(cifogParam, snap, thisRedshiftList->redshifts[snap], deltaRedshift, grid, cifogNion, cifogNionHeI, cifogNionHeII, integralTable, photIonBgList, thisRank);
      
      /* update current photoionization & XHII field */
      /* WATCH OUT: The following two lines are order sensitive! */
      update_field(&photHI, grid, thisRedshiftList->redshifts[snap], simParam->XHII_threshold, photHIName, memory_intensive);
      update_field(&zreion, grid, thisRedshiftList->redshifts[snap], simParam->XHII_threshold, zreionName, memory_intensive);
      
      /* get photoionization and ionization values at galaxies' positions */
      if(simParam->radfeedback == 1 && memory_intensive == 0)
        map_grid_to_gal(numGtrees, &theseGtrees, minSnap, maxSnap, simParam, domain, grid);

      /* save reionization redshift grid in a file */
      write_zreion_grids_to_file(simParam, domain, snap, grid, zreion, photHI);
      
      /* initialise ionizing emissivity list */
      if(cifogNion != NULL) fftw_free(cifogNion);
      deallocate_nion(thisNionList); 
    }
    
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
        
    /* save the number of all galaxies in an array */
    if(simParam->horizontalOutput == 1)
    {
      int printNumGalSnap = numGalSnap[snap];
#ifdef MPI
      MPI_Allreduce(&(numGalSnap[snap]), &printNumGalSnap, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
      if(thisRank == 0) printf("snap %d: numGalSnap = %d\n", snap, printNumGalSnap);
    }
  }

  /* Write number of galaxies contained in each snap into a file */
  if(simParam->horizontalOutput == 1)
  {
    write_num_galaxies_to_file(simParam, numGalSnap, endSnap, thisRank);
  }
  
  /* Write trees into #PROC files */
  if(simParam->verticalOutput == 1)
  {
    if(simParam->endSnap!=OUTPUTSNAPNUMBER+1)
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
  long int tmpNumHalos = numHalos;
  MPI_Allreduce(&tmpNumHalos, &numHalos, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(thisRank == 0) printf("Total number of halos/galaxies = %ld\n", numHalos);
#endif

  /*---------------------------------------------------------------------*/
  /* DEALLOCATION */
  /*---------------------------------------------------------------------*/

  /* end CIFOG */
  if(simParam->reion == 1 && simParam->reion_model == 1)
  {
      cifog_deallocate(cifogParam, grid, integralTable, photIonBgList, thisRank);
      if(photHI != NULL) fftw_free(photHI);
      if(zreion != NULL) fftw_free(zreion);
  }

  /* deallocate all lists, trees and other structs */
  deallocate_redshiftlist(thisRedshiftList);
  free(numGalSnap);
  deallocate_gtreeList(theseGtrees, numGtrees);
  deallocate_domain(domain);
}

/*---------------------------------------------------------------------*/
/* This function evolves all galaxies in the trees by (maxSnap-minSnap)*/
/* time steps                                                          */
/*---------------------------------------------------------------------*/
void evolve(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, dconfObj_t simParam, domain_t *thisDomain, nion_t *thisNionList, fftw_complex *photHI, fftw_complex *zreion, int outSnap, outgalsnap_t **outGalList, int *numGalInSnapList)
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
            /* get grid properties for memory intensive option (zreion and photHI grids exist) */
            if(zreion != NULL && simParam->reion == 1)
            {
              thisGal->photHI_bg = get_gridProperty_gal(thisGal, thisDomain->nbins, inv_boxsize, photHI);

              if(thisGal->zreion <= 0.) 
              {
                thisGal->zreion = get_gridProperty_gal(thisGal, thisDomain->nbins, inv_boxsize, zreion);
              }
            }
            
            /* evolve galaxy */
            evolve_gal(thisGal, thisGtree, simParam);
            
            if(thisGal->zreion > 0.) assert(thisGal->zreion >= 1./thisGal->scalefactor - 1.);

            /* store ionizing properties for reionization */
            if(simParam->reion == 1) copy_gal_to_nion(thisGal, simParam, thisDomain, 1.1, thisNionList);
            
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
                for(int outGal=numGalInOutSnap; outGal<MAXGAL + count * MAXGAL_INC; outGal++)
                  initOutGalSnapWithoutAllocation(&((*outGalList)[outGal]));
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

/*---------------------------------------------------------------------*/
/* This function sets the pointer within each tree to the first galaxy */
/* at snapshot <startSnap>                                             */
/*---------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------*/
/* This function writes the grids containing the reionization redshift */
/* and the photoionization rate at that redshift to files              */
/*---------------------------------------------------------------------*/
void write_zreion_grids_to_file(dconfObj_t simParam, domain_t *thisDomain, int snap, grid_t *thisGrid, fftw_complex *zreion, fftw_complex *photHI)
{
  int memory_intensive = simParam->memoryIntensive;
  char zreionFile[MAXLENGTH];
  for(int i=0; i<MAXLENGTH; i++) zreionFile[i] = '\0';
  sprintf(zreionFile, "%s/zreion_%02d", simParam->outFileName, snap);
  if(memory_intensive == 1)
  {
    write_grid_to_file_double_test(zreion, thisGrid->nbins, thisGrid->nbins, 0, zreionFile, thisDomain->thisRank, memory_intensive);
  }
  else
  {
    write_grid_to_file_double_test(thisGrid->zreion, thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, zreionFile, thisDomain->thisRank, memory_intensive);
  }
  
  /* save photoionization grid in a file */
  char photHIFile[MAXLENGTH];
  for(int i=0; i<MAXLENGTH; i++) photHIFile[i] = '\0';
  sprintf(photHIFile, "%s/photHI_zreion_%02d", simParam->outFileName, snap);
  if(memory_intensive == 1)
  {
    write_grid_to_file_double_test(photHI, thisGrid->nbins, thisGrid->nbins, 0, photHIFile, thisDomain->thisRank, memory_intensive);
  }
  else
  {
    write_grid_to_file_double_test(thisGrid->photHI_zreion, thisGrid->nbins, thisGrid->local_n0, thisGrid->local_0_start, photHIFile, thisDomain->thisRank, memory_intensive);
  }
}
