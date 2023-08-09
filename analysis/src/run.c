#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "dconfObj.h"
#include "outgal.h"
#include "read_trees.h"
#include "zas_lists.h"

#include "build_index_tree_walking.h"
#include "derive_properties.h"
#include "statistics.h"
#include "selection.h"
#include "selection_median.h"
#include "selection_evolution.h"
#include "numdens_analytic.h"
#include "selection_analytic.h"
#include "selection_evolution_analytic.h"
#include "selection_galaxypairs.h"
#include "tree_correct.h"

#include "SNfeedback.h"
#include "ion_emissivity.h"

#include "output_lists.h"

void analysis(dconfObj_t simParam, int thisRank, int size)
{
  if(directory_exist(simParam->outputDir) == 0)
  {
    mkdir(simParam->outputDir, 0744);
  }
  
  /* -----------------------------------------------------------------------*/
  /* LOADING TREES */
  /* -----------------------------------------------------------------------*/
  outgtree_t **theseTrees = NULL;
  int32_t numTrees = read_trees_in_all_files(simParam->inputFile, simParam->inputType, simParam->numFiles, &theseTrees, thisRank, size);
    
  /* -----------------------------------------------------------------------*/
  /* getting snapnumbers, scalefactors, redshifts & times */
  /* -----------------------------------------------------------------------*/
  simParam->numSnaps = get_num_snapshots(theseTrees);
  simParam->scalefactors = get_scalefactors(theseTrees, numTrees); 
  simParam->redshifts = get_redshifts(simParam->scalefactors, simParam->numSnaps);
  simParam->times = get_times_from_redshifts(simParam->redshifts, simParam->numSnaps, simParam->hubble_h, simParam->omega_m, simParam->omega_l);
  simParam->SNenergy = get_SNenergy_delayed(simParam->numSnaps, simParam->times);
  simParam->corrFactor_nion = get_corrFactor_nion(simParam, simParam->numSnaps);

  correct_localDescID_in_trees(theseTrees, numTrees, simParam->numSnaps-1);
  
  /* -----------------------------------------------------------------------*/
  /* BUILD INDEX FOR TREE WALKING */
  /* -----------------------------------------------------------------------*/
  int **index = malloc(numTrees * sizeof(int*));
  for(int tree=0; tree<numTrees; tree++)
  {
    index[tree] = malloc(theseTrees[tree]->numGal*sizeof(int));
    for(int gal=0; gal<theseTrees[tree]->numGal; gal++)
    {
      index[tree][gal] = -1;
    }
  }
                  
  int ***listEquals = malloc(numTrees * sizeof(int**));
  int *sizeListEquals = malloc(numTrees * sizeof(int));
  int **listMerged = malloc(numTrees * sizeof(int*));
  for(int tree=0; tree<numTrees; tree++)
  {
    sizeListEquals[tree] = 0;
    listEquals[tree] = NULL;
    listMerged[tree] = init_toBeMerged();
  }
  
  int prevSnap = -1;
  int currSnap = 0;
  int endSnap = simParam->numSnaps-1;
  if(simParam->numOutputRedshifts > 0)
    endSnap = find_snapshoft_from_redshift(simParam->outputRedshifts[simParam->numOutputRedshifts-1], simParam->numSnaps, simParam->redshifts);
  else if(simParam->trackEvolution != 1)
  {
    fprintf(stderr, "No output redshift is given. Aborting...\n");
    exit(EXIT_FAILURE);
  }
        
  if(simParam->numOutputRedshifts > 0 && find_snapshoft_from_redshift(simParam->outputRedshifts[0], simParam->numSnaps, simParam->redshifts) <= endSnap)
  {
    if(thisRank == 0) printf("endSnap = %d corresponds to z = %f\n", endSnap, simParam->redshifts[endSnap]);
    int counter = 0;
    
    for(int snap=0; snap<=endSnap; snap++)
    {   
      currSnap = snap;
      get_listEquals(theseTrees, numTrees, index, listEquals, sizeListEquals, listMerged, prevSnap, currSnap);
      prevSnap = currSnap;
      
      if(snap == find_snapshoft_from_redshift(simParam->outputRedshifts[counter], simParam->numSnaps, simParam->redshifts))
      {
        if(thisRank == 0) printf("outputing at z = %f or snap = %d\n", simParam->outputRedshifts[counter], find_snapshoft_from_redshift(simParam->outputRedshifts[counter], simParam->numSnaps, simParam->redshifts));
        
  #ifdef MPI
        MPI_Barrier(MPI_COMM_WORLD);
  #endif

        /* -----------------------------------------------------------------------*/
        /* DERIVING & ANALYSING PROPERTIES */
        /* -----------------------------------------------------------------------*/
        
        for(int i=0; i<simParam->num_2D_history; i++)
        {
          char *filename = create_filename_2D_histogram(simParam->outputDir, simParam->property_2D_history[i], simParam->binProperty1_2D_history[i], simParam->binProperty2_2D_history[i], simParam->outputRedshifts[counter], simParam->smoothingScale, simParam->MvirThreshold);
          
          get_2D_histogram_history(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->property_2D_history[i], simParam->binProperty1_2D_history[i], simParam->binProperty2_2D_history[i], simParam->binsInLog1_2D_history[i], simParam->binsInLog2_2D_history[i], simParam->binsPerMag1_2D_history[i], simParam->binsPerMag2_2D_history[i], 0, 0, filename, thisRank);
          free(filename);
        }
        
        for(int i=0; i<simParam->num_1D_history; i++)
        {
          char *filename = create_filename_1D_histogram(simParam->outputDir, simParam->property_1D_history[i], simParam->binProperty_1D_history[i], simParam->outputRedshifts[counter], simParam->smoothingScale, simParam->MvirThreshold);
          
          get_1D_histogram_history(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->property_1D_history[i], simParam->binProperty_1D_history[i], simParam->binsInLog_1D_history[i], simParam->binsPerMag_1D_history[i], 0, filename, thisRank);
          free(filename);
        }
        
        for(int i=0; i<simParam->num_2D_history_median; i++)
        {
          char *filename = create_filename_2D_histogram_history_median(simParam->outputDir, simParam->property_2D_history_median[i], simParam->binProperty1_2D_history_median[i], simParam->binProperty2_2D_history_median[i], simParam->outputRedshifts[counter], simParam->smoothingScale, simParam->MvirThreshold);
          
          get_2D_histogram_history_median(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->property_2D_history_median[i], simParam->binProperty1_2D_history_median[i], simParam->binProperty2_2D_history_median[i], simParam->binsInLog1_2D_history_median[i], simParam->binsInLog2_2D_history_median[i], simParam->binsPerMag1_2D_history_median[i], simParam->binsPerMag2_2D_history_median[i], 0, 0, filename, thisRank, size);
          free(filename);
        }
        
        for(int i=0; i<simParam->num_3D_value; i++)
        {
          char *filename = create_filename_3D_histogram_value(simParam->outputDir, simParam->property_3D_value[i], simParam->binProperty1_3D_value[i], simParam->binProperty2_3D_value[i], simParam->binProperty3_3D_value[i], simParam->outputRedshifts[counter], simParam->smoothingScale, simParam->MvirThreshold);
          
          get_3D_histogram_value(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->property_3D_value[i], simParam->binProperty1_3D_value[i], simParam->binProperty2_3D_value[i], simParam->binProperty3_3D_value[i], simParam->binsInLog1_3D_value[i], simParam->binsInLog2_3D_value[i], simParam->binsInLog3_3D_value[i], simParam->binsPerMag1_3D_value[i], simParam->binsPerMag2_3D_value[i], simParam->binsPerMag3_3D_value[i], 0, 0, 0, (float)(simParam->binProperty1_3D_mapLowLimit[i]), (float)(simParam->binProperty1_3D_mapUpLimit[i]), filename, thisRank);
          free(filename);
        }
        
        for(int i=0; i<simParam->num_2D_value; i++)
        {
          char *filename = create_filename_2D_histogram_value(simParam->outputDir, simParam->property_2D_value[i], simParam->binProperty1_2D_value[i], simParam->binProperty2_2D_value[i], simParam->outputRedshifts[counter], simParam->smoothingScale, simParam->MvirThreshold);
          
          get_2D_histogram_value(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->property_2D_value[i], simParam->binProperty1_2D_value[i], simParam->binProperty2_2D_value[i], simParam->binsInLog1_2D_value[i], simParam->binsInLog2_2D_value[i], simParam->binsPerMag1_2D_value[i], simParam->binsPerMag2_2D_value[i], 0, 0, (float)(simParam->binProperty1_2D_mapLowLimit[i]), (float)(simParam->binProperty1_2D_mapUpLimit[i]), filename, thisRank);
          free(filename);
        }

        for(int i=0; i<simParam->num_1D_value; i++)
        {
          char *filename = create_filename_1D_histogram_value(simParam->outputDir, simParam->property_1D_value[i], simParam->binProperty_1D_value[i], simParam->outputRedshifts[counter], simParam->smoothingScale, simParam->MvirThreshold);
          
          get_1D_histogram_value(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->property_1D_value[i], simParam->binProperty_1D_value[i], simParam->binsInLog_1D_value[i], simParam->binsPerMag_1D_value[i], (float)(simParam->binProperty1_1D_mapLowLimit[i]), 0, (float)(simParam->binProperty1_1D_mapUpLimit[i]), filename, thisRank);
          free(filename);
        }
        
        for(int i=0; i<simParam->num_2D_median; i++)
        {
          char *filename = create_filename_2D_histogram_median(simParam->outputDir, simParam->property_2D_median[i], simParam->binProperty1_2D_median[i], simParam->binProperty2_2D_median[i], simParam->outputRedshifts[counter], simParam->smoothingScale, simParam->MvirThreshold);
          
          get_2D_histogram_median(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->property_2D_median[i], simParam->binProperty1_2D_median[i], simParam->binProperty2_2D_median[i], simParam->binsInLog1_2D_median[i], simParam->binsInLog2_2D_median[i], simParam->binsPerMag1_2D_median[i], simParam->binsPerMag2_2D_median[i], 0, 0, (float)(simParam->binProperty1_2D_median_mapLowLimit[i]), (float)(simParam->binProperty1_2D_median_mapUpLimit[i]), filename, thisRank, size);
          free(filename);
        }
    
        for(int i=0; i<simParam->num_2D; i++)
        {
          char *filename = create_filename_2D_histogram_numDens(simParam->outputDir, simParam->binProperty1_2D[i], simParam->binProperty2_2D[i], simParam->outputRedshifts[counter], simParam->smoothingScale, simParam->MvirThreshold);
          
          if(simParam->analytic == 1)
          {
            char *hmfFilename = create_filename_hmf(simParam->hmfFilename, simParam->mergertreeEndRedshift);
            get_2D_histogram_analytic(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->binProperty1_2D[i], simParam->binProperty2_2D[i], simParam->binsInLog1_2D[i], simParam->binsInLog2_2D[i], simParam->binsPerMag1_2D[i], simParam->binsPerMag2_2D[i], 0, 0, filename, hmfFilename, thisRank);
            free(hmfFilename);
          }
          else
            get_2D_histogram(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->binProperty1_2D[i], simParam->binProperty2_2D[i], simParam->binsInLog1_2D[i], simParam->binsInLog2_2D[i], simParam->binsPerMag1_2D[i], simParam->binsPerMag2_2D[i], 0, 0, filename, thisRank);
          free(filename);
        }
        
        for(int i=0; i<simParam->num_1D; i++)
        {
          char *filename = create_filename_1D_histogram_numDens(simParam->outputDir, simParam->binProperty_1D[i], simParam->outputRedshifts[counter], simParam->smoothingScale, simParam->MvirThreshold);
          
          if(simParam->analytic == 1)
          {
            char *hmfFilename = create_filename_hmf(simParam->hmfFilename, simParam->mergertreeEndRedshift);
            get_1D_histogram_analytic(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->binProperty_1D[i], simParam->binsInLog_1D[i], simParam->binsPerMag_1D[i], 0, simParam->cumulative[i], filename, hmfFilename, thisRank);
            free(hmfFilename);
          }
          else
            get_1D_histogram(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->binProperty_1D[i], simParam->binsInLog_1D[i], simParam->binsPerMag_1D[i], 0, simParam->cumulative[i], filename, thisRank);
          free(filename);
        }
        
        for(int i=0; i<simParam->num_galaxyPairs; i++)
        {
          char *filename = create_filename_galaxypairs(simParam->outputDir, simParam->selectionProperty1[i], simParam->selectionProperty2[i], simParam->minSelectionProperty1[i], simParam->maxSelectionProperty1[i], simParam->minSelectionProperty2[i], simParam->maxSelectionProperty2[i], simParam->outputRedshifts[counter]);

          select_galaxypairs(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, simParam->selectionProperty1[i], simParam->selectionProperty2[i], simParam->propertyWithHistory_galaxyPairs[i], simParam->minSelectionProperty1[i], simParam->maxSelectionProperty1[i], simParam->minSelectionProperty2[i], simParam->maxSelectionProperty2[i], simParam->maxDistanceInComMpc[i], filename, thisRank, size);
          free(filename);
        }
        
        counter++;
        if(counter == simParam->numOutputRedshifts)
          break;
      }
    }
  }
  
  /* -----------------------------------------------------------------------*/
  /* DEALLOCATION */
  /* -----------------------------------------------------------------------*/
  
  for(int tree=0; tree<numTrees; tree++)
  {
    free(index[tree]);
    for(int gal=0; gal<sizeListEquals[tree]; gal++)
      free(listEquals[tree][gal]);
    free(listEquals[tree]);
    free(listMerged[tree]);
  }
  free(index);
  free(listEquals);
  free(sizeListEquals);
  free(listMerged);
  
  
  if(simParam->trackEvolution == 1)
  {
    for(int i=0; i<simParam->num_1D_evolution; i++)
    {
      /* -----------------------------------------------------------------------*/
      /* REALLOCATION */
      /* -----------------------------------------------------------------------*/
      
      index = malloc(numTrees * sizeof(int*));
      for(int tree=0; tree<numTrees; tree++)
      {
        index[tree] = malloc(theseTrees[tree]->numGal*sizeof(int));
        for(int gal=0; gal<theseTrees[tree]->numGal; gal++)
          index[tree][gal] = -1;
      }
                  
      listEquals = malloc(numTrees * sizeof(int**));
      sizeListEquals = malloc(numTrees * sizeof(int));
      listMerged = malloc(numTrees * sizeof(int*));
      for(int tree=0; tree<numTrees; tree++)
      {
        sizeListEquals[tree] = 0;
        listEquals[tree] = NULL;
        listMerged[tree] = init_toBeMerged();
      }
      
      /* -----------------------------------------------------------------------*/
      /* TRACKING EVOLUTION */
      /* -----------------------------------------------------------------------*/
      
      if(simParam->analytic == 1)
      {
        get_1D_histogram_evolution_analytic(simParam, theseTrees, numTrees, index, listEquals, sizeListEquals, listMerged, i, thisRank);
      }
      else
        get_1D_histogram_evolution(simParam, theseTrees, numTrees, index, listEquals, sizeListEquals, listMerged, i, thisRank);
      
      /* -----------------------------------------------------------------------*/
      /* DEALLOCATION */
      /* -----------------------------------------------------------------------*/
      
      for(int tree=0; tree<numTrees; tree++)
      {
        free(index[tree]);
        for(int gal=0; gal<sizeListEquals[tree]; gal++)
          free(listEquals[tree][gal]);
        free(listEquals[tree]);
        free(listMerged[tree]);
      }
      free(index);
      free(listEquals);
      free(sizeListEquals);
      free(listMerged);
    }
  }
  
  if(simParam->outputLists == 1 && simParam->analytic == 0)
  {
    /* -----------------------------------------------------------------------*/
    /* REALLOCATION */
    /* -----------------------------------------------------------------------*/
    
    index = malloc(numTrees * sizeof(int*));
    for(int tree=0; tree<numTrees; tree++)
    {
      index[tree] = malloc(theseTrees[tree]->numGal*sizeof(int));
      for(int gal=0; gal<theseTrees[tree]->numGal; gal++)
        index[tree][gal] = -1;
    }
                
    listEquals = malloc(numTrees * sizeof(int**));
    sizeListEquals = malloc(numTrees * sizeof(int));
    listMerged = malloc(numTrees * sizeof(int*));
    for(int tree=0; tree<numTrees; tree++)
    {
      sizeListEquals[tree] = 0;
      listEquals[tree] = NULL;
      listMerged[tree] = init_toBeMerged();
    }
    
    /* -----------------------------------------------------------------------*/
    /* WRITE ALL GALAXIES */
    /* -----------------------------------------------------------------------*/
    
    write_galaxies_to_textfile_all_snaps(simParam, theseTrees, numTrees, index, listEquals, sizeListEquals, listMerged, thisRank);
    
    /* -----------------------------------------------------------------------*/
    /* DEALLOCATION */
    /* -----------------------------------------------------------------------*/
    
    for(int tree=0; tree<numTrees; tree++)
    {
      free(index[tree]);
      for(int gal=0; gal<sizeListEquals[tree]; gal++)
        free(listEquals[tree][gal]);
      free(listEquals[tree]);
      free(listMerged[tree]);
    }
    free(index);
    free(listEquals);
    free(sizeListEquals);
    free(listMerged);
  }
  
  deallocate_outgtreeList(theseTrees, numTrees);
}
