#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <complex.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "../utils.h"

#include "phys_const.h"
#include "confObj.h"
#include "grid.h"
#include "photion_background.h"
#include "fraction_q.h"
#include "filtering.h"
#include "self_shielding.h"

#include "density_distribution.h"
#include "recombination.h"
#include "mean_free_path.h"

#include "input_redshifts.h"
#include "input_grid.h"
#include "input_nion.h"
#include "checks.h"
#include "restart.h"

#include "cifog_delphi.h"
 
//-------------------------------------------------------------------------------
// reading input files and prepare grid
//-------------------------------------------------------------------------------
    
int cifog_init(char *iniFile, confObj_t *simParam, grid_t **grid,  integral_table_t **integralTable, photIonlist_t **photIonBgList, const int myRank)
{
    if(myRank==0) printf("\n********************\nINITIALIZING CIFOG\n********************\n");
    
#ifdef MPI
    fftw_mpi_init();
#endif
    
    //read paramter file
    *simParam = readConfObj(iniFile);
    
    //verify that helium runs contain helium!
    if((*simParam)->solve_He == 1)
    {
        if((*simParam)->Y <= 0.)
        {
            fprintf(stderr, "If you include helium its mass fraction Y should be larger than zero!\n");
            exit(EXIT_FAILURE);
        }
    }
    
    check_output_directories_exist((*simParam));
    
    //read files (allocate grid)
    *grid = initGrid();
    if(myRank==0) printf("\n++++\nInitialising grids... ");
    read_files_to_grid(*grid, (*simParam));
    if(myRank==0) printf("done\n+++\n");
    
    //read photoionization background values 
    if((*simParam)->photHI_model == 11)
    {
        if(myRank==0) printf("\n++++\nReading photoionization background rates... ");
        *photIonBgList = read_photIonlist((*simParam)->photHI_bg_file);
        if(myRank==0) printf("done\n+++\n");
    }
    
    if((*simParam)->calc_recomb == 2)
    {
        //read table for recombinations
        if(myRank==0) printf("\n++++\nRead table for recombinations... ");
        *integralTable = initIntegralTable((*simParam)->zmin, (*simParam)->zmax, (*simParam)->dz, (*simParam)->fmin, (*simParam)->fmax, (*simParam)->df, (*simParam)->dcellmin, (*simParam)->dcellmax, (*simParam)->ddcell);
        if(myRank==0) printf("done\n+++\n");
    }

    return EXIT_SUCCESS;
}


int cifog(confObj_t simParam, int snap, double redshift, double deltaRedshift, grid_t *grid, fftw_complex *nion, fftw_complex *nion_HeI, fftw_complex *nion_HeII, const integral_table_t *integralTable, photIonlist_t *photIonBgList, const int myRank)
{
  /* FIX REDSHIFT & EVOLUTION TIME */
  simParam->redshift_prev_snap = redshift;
  simParam->redshift = redshift - deltaRedshift;
  
  int cycle_offset = simParam->snapshot_start;
  int cycle  = snap - cycle_offset;

  if(myRank==0) printf("snap = %d\t cycle_offset  = %d\t cycle = %d\n", snap, cycle_offset, cycle);
  
  cifog_step(simParam, grid, nion, nion_HeI, nion_HeII, integralTable, photIonBgList, cycle, cycle_offset, snap, myRank);
  
  return EXIT_SUCCESS;
}

int cifog_step(confObj_t simParam, grid_t *grid, fftw_complex *nion, fftw_complex *nion_HeI, fftw_complex *nion_HeII, const integral_table_t *integralTable, photIonlist_t *photIonBgList, const int cycle, const int cycle_offset, int snap, const int myRank)
{
    const double f = simParam->f;
    char photHIFile[MAXLENGTH], XionFile[MAXLENGTH];
    
    if(myRank==0)
    {
        printf("\n******************\nSNAP %d\t CYCLE %d\n******************\n", snap, cycle);
    }
    
    if(myRank==0) printf("\n++++\nReading sources/nion file for snap = %d... ", snap);
    read_update_nion(simParam, nion, grid);
    if(myRank==0) printf("done\n+++\n");
    
//     //write nion field to file
//     char nionFile[MAXLENGTH];
//     for(int i=0; i<MAXLENGTH; i++) nionFile[i] = '\0';
//     sprintf(nionFile, "%s_%02d", "/net/eos/data/users/hutter/test_delphi/nion.out", cycle + cycle_offset);
//     if(myRank==0) printf("\n++++\nwriting nion field to file %s ... ", nionFile);
//     save_to_file(grid->nion, grid, nionFile);
//     if(myRank==0) printf("done\n+++\n");
    
    if(simParam->solve_He == 1)
    {
        if(myRank==0) printf("\n++++\nReading sources/nion file for snap = %d... ", snap);
        read_update_nion_HeI(simParam, nion_HeI, grid);
        if(myRank==0) printf("done\n+++\n");
        
        if(myRank==0) printf("\n++++\nReading sources/nion file for snap = %d... ", snap);
        read_update_nion_HeII(simParam, nion_HeII, grid);
        if(myRank==0) printf("done\n+++\n");
    }
    
    if(myRank==0) printf("\n++++\nReading igm density file for snap = %d... ", snap);
    read_update_igm_density(simParam, grid, snap);
    if(myRank==0) printf("done\n+++\n");
  
    if(myRank==0) printf("\n++++\nReading igm clump file for snap = %d... ", snap);
    read_update_igm_clump(simParam, grid, snap);
    if(myRank==0) printf("done\n+++\n");
    
    //------------------------------------------------------------------------------
    // compute web model
    //------------------------------------------------------------------------------
    
    if(simParam->use_web_model == 1)
    {
        if(cycle == 0)
        {
            /* ----------------------------------------- */
            /* photoionization rate is homogeneous       */
            /* ----------------------------------------- */
            if(simParam->photHI_model == 0)
            {
                //set photoionization rate on grid to background value
                if(myRank==0) printf("\n++++\nPHOTHI_CONST: setting photoionization rate to background value...");
                set_value_to_photoionization_field(grid, simParam);
                if(myRank==0) printf("\n photHI_bg = %e s^-1\n", simParam->photHI_bg);
                if(myRank==0) printf("done\n+++\n");
            }
            
            /* ----------------------------------------------------------------------------------------- */
            /* photoionization rate is given by mean mfp and exp(-r/mfp)/r^2 around sources distribution */
            /* ----------------------------------------------------------------------------------------- */
            else if(simParam->photHI_model == 11)
            {
                if(myRank==0) printf("\n++++\nPHOTHI_GIVEN: computing photoionization rate...");
                //set photoionization value according to the given list
                set_value_to_photHI_bg(simParam, get_photHI_from_redshift(photIonBgList, simParam->redshift));
                
                if(simParam->calc_mfp == 1)
                {
                    simParam->mfp = f*simParam->box_size/(simParam->h * (1.+simParam->redshift))/grid->nbins;
                    if(myRank==0) printf("\n mean free path at z = %e is %e Mpc\n", simParam->redshift, simParam->mfp);
                }
                
                //compute spatial photoionization rate according to source distribution
                compute_photHI(grid, simParam, 1);
                if(myRank==0) printf("done\n+++\n");
            }
            
            /* ----------------------------------------------------------------------------------------- */
            /* photoionization rate is given by mean mfp and exp(-r/mfp)/r^2 around sources distribution */
            /* ----------------------------------------------------------------------------------------- */
            else if(simParam->photHI_model == 1)
            {                
                if(myRank==0) printf("\n++++\nPHOTHI_FLUX: computing photoionization rate...");

                if(simParam->calc_mfp == 1)
                {
                    simParam->mfp = f*simParam->box_size/(simParam->h * (1.+simParam->redshift))/grid->nbins;
                    if(myRank==0) printf("\n mean free path at z = %e is %e Mpc\n", simParam->redshift, simParam->mfp);
                }
                
                //compute spatial photoionization rate according to source distribution
                compute_photHI(grid, simParam, 0);
                if(myRank==0) printf("done\n+++\n");
            }
            
            /* ------------------------------------------------------- */
            /* photoionization rate is given by mfp of ionized regions */
            /* ------------------------------------------------------- */
            else if(simParam->photHI_model == 2)
            {                
                if(myRank==0) printf("\n++++\nPHOTHI_MFP: set photoionization rate according to ionized regions... ");
                if(cycle != 0){
                    compute_photHI_ionizedRegions(grid, simParam);
                }else{
                    set_value_to_photoionization_field(grid,simParam);
                }
                if(myRank==0) printf("done\n+++\n");
            }
            
            else
            {
                if(myRank==0)
                {
                    printf("\n+++\nNo supported photoionization rate model. Photoionization model is required for the web model. Abborting...\n");
                    exit(EXIT_FAILURE);
                }
            }
        }
        
        /* ------------------------------------------------------- */
        /* compute HI fraction (web model)                         */
        /* ------------------------------------------------------- */
        if(myRank==0) printf("\n++++\nApply web model... ");
        compute_web_ionfraction(grid, simParam);
        if(myRank==0) printf("done\n+++\n");
        
        
        /* ------------------------------------------------------- */
        /* DISABLED: compute local mfp in each cell                */
        /* ------------------------------------------------------- */
        if(simParam->calc_mfp == -1)
        {
            //compute mean free paths
            if(myRank==0) printf("\n++++\nCompute mean free paths... ");
            compute_web_mfp(grid, simParam);
            if(myRank==0) printf("done\n+++\n");
        }
    }
    
    /* ---------------------------------------------------------------- */
    /* compute recombinations based on local photion rate & HI fraction */
    /* ---------------------------------------------------------------- */
    if(simParam->calc_recomb == 1 && simParam->const_recomb == 0)
    {
        //compute number of recombinations (HII, HeII & HeIII)
        if(myRank==0) printf("\n++++\nCompute number of recombinations... ");
        compute_number_recombinations(grid, simParam);
        if(myRank==0) printf("done\n+++\n");
    }
    else if(simParam->calc_recomb == 2 && simParam->const_recomb == 0)
    {
        //compute number of recombinations (HII, HeII & HeIII)
        if(myRank==0) printf("\n++++\nCompute number of recombinations... ");
        compute_number_recombinations_M2000(grid, simParam, simParam->recomb_table, integralTable);
        if(myRank==0) printf("done\n+++\n");
    }
    
    //------------------------------------------------------------------------------
    // compute number of recombinations for a constant rate
    //------------------------------------------------------------------------------
    
    if(simParam->const_recomb == 1)
    {
        //compute number of recombinations (HII, HeII & HeIII)
        if(myRank==0) printf("\n++++\nCompute number of recombinations... ");
        compute_number_recombinations_const(grid, simParam);
        if(myRank==0) printf("done\n+++\n");
    }

    //--------------------------------------------------------------------------------
    // apply tophat filter
    //--------------------------------------------------------------------------------
    
    //compute fraction Q
    if(myRank==0) printf("\n++++\nHII: computing relation between number of ionizing photons and absorptions... ");
    compute_cum_values(grid, simParam, 0, myRank);
    if(myRank==0) printf("done\n+++\n");
    
    //apply filtering
    if(myRank==0) printf("\n++++\nHII: apply tophat filter routine for ionization field... ");
    compute_ionization_field(simParam, grid, photIonBgList, 0, myRank);
    if(myRank==0) printf("done\n+++\n");
    
    //write ionization field to file
    for(int i=0; i<MAXLENGTH; i++) XionFile[i] = '\0';
    sprintf(XionFile, "%s_%02d", simParam->out_XHII_file, cycle + cycle_offset);
    if(myRank==0) printf("\n++++\nwriting HI ionization field to file %s ... ", XionFile);
    save_to_file(grid->XHII, grid, XionFile);
    if(myRank==0) printf("done\n+++\n");
    
    if(simParam->solve_He == 1)
    {
        //compute fraction Q
        if(myRank==0) printf("\n++++\nHeII/HeIII: computing relation between number of ionizing photons and absorptions... ");
        compute_cum_values(grid, simParam, 1, myRank);
        compute_cum_values(grid, simParam, 2, myRank);
        if(myRank==0) printf("done\n+++\n");
    
        //apply filtering
        if(myRank==0) printf("\n++++\nHeII/HeIII: apply tophat filter routine for ionization field... ");
        compute_ionization_field(simParam, grid, photIonBgList, 1, myRank);
        compute_ionization_field(simParam, grid, photIonBgList, 2, myRank);
        if(myRank==0) printf("done\n+++\n");
        
        //write ionization field to file
        for(int i=0; i<MAXLENGTH; i++) XionFile[i] = '\0';
        sprintf(XionFile, "%s_%02d", simParam->out_XHeII_file, cycle + cycle_offset);
        if(myRank==0) printf("\n++++\nWriting HeI ionization field to file %s ... ", XionFile);
        save_to_file(grid->XHeII, grid, XionFile);
        if(myRank==0) printf("done\n+++\n");
    
        //write ionization field to file
        for(int i=0; i<MAXLENGTH; i++) XionFile[i] = '\0';
        sprintf(XionFile, "%s_%02d", simParam->out_XHeIII_file, cycle + cycle_offset);
        if(myRank==0) printf("\n++++\nWriting HeII ionization field to file %s ... ", XionFile);
        save_to_file(grid->XHeIII, grid, XionFile);
        if(myRank==0) printf("done\n+++\n");
    }
    
    if(simParam->use_web_model == 1 && simParam->write_photHI_file == 1)
    {
        //write photoionization rate field to file
        for(int i=0; i<MAXLENGTH; i++) photHIFile[i] = '\0';
        sprintf(photHIFile, "%s_%02d", simParam->out_photHI_file, cycle + cycle_offset);
        if(myRank==0) printf("\n++++\nWriting HI photoionization field to file... ");
        save_to_file(grid->photHI, grid, photHIFile);
        if(myRank==0) printf("done\n+++\n");
    }
            
    return EXIT_SUCCESS;
}

int cifog_deallocate(confObj_t simParam, grid_t *grid, integral_table_t *integralTable, photIonlist_t *photIonBgList, const int myRank)
{
    //--------------------------------------------------------------------------------
    // deallocating grids
    //--------------------------------------------------------------------------------
  
    if(myRank==0) printf("\n********************\nFINALIZING CIFOG\n********************\n");

    if(simParam->calc_recomb == 2)
    {
        //read table for recombinations
        if(myRank==0) printf("\n++++\nDeallocating table for recominsations... ");
        free(integralTable);
        integralTable  = NULL;
        if(myRank==0) printf("done\n+++\n");
    }

    if(myRank==0) printf("\n++++\nDeallocating background photionization rate list... ");
    deallocate_photIonlist(photIonBgList);
    photIonBgList = NULL;
    if(myRank==0) printf("done\n+++\n");

    //deallocate grid
    if(myRank==0) printf("\n++++\nDeallocating grid ...");
    deallocate_grid(grid, simParam);
    grid = NULL;
    if(myRank==0) printf("done\n+++\n");
    
    confObj_del(&simParam);
    
#ifdef MPI
    fftw_mpi_cleanup();
#endif
    
    return EXIT_SUCCESS;
}