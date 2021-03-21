/*
 *  dconfObj.c
 *  uvff
 *
 *  Created by 
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

/*--- Includes ----------------------------------------------------------*/
#include "dconfObj.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>
#include "xmem.h"
#include "utils.h"

/*--- Defines for the Ini structure -------------------------------------*/


/*--- Prototypes of local functions -------------------------------------*/


/*--- Implementations of exported functions ------------------------------*/
extern dconfObj_t


dconfObj_new(parse_ini_t ini)
{
    dconfObj_t config;
    assert(ini != NULL);
    
    config = xmalloc(sizeof(struct dconfObj_struct));
    
    config->inputType = 1;
    config->outputRedshifts = NULL;
    
    config->analytic = 0;
    config->mergertreeBinwidth = 0.;
    config->hmfFilename = NULL;
    
    config->property_2D_history = NULL;
    config->binProperty1_2D_history = NULL;
    config->binProperty2_2D_history = NULL;
    config->binsInLog1_2D_history = NULL;
    config->binsInLog2_2D_history = NULL;
    config->binsPerMag1_2D_history = NULL;
    config->binsPerMag2_2D_history = NULL;

    config->property_1D_history = NULL;
    config->binProperty_1D_history = NULL;
    config->binsInLog_1D_history = NULL;
    config->binsPerMag_1D_history = NULL;
    
    config->property_2D_history_median = NULL;
    config->binProperty1_2D_history_median = NULL;
    config->binProperty2_2D_history_median = NULL;
    config->binsInLog1_2D_history_median = NULL;
    config->binsInLog2_2D_history_median = NULL;
    config->binsPerMag1_2D_history_median = NULL;
    config->binsPerMag2_2D_history_median = NULL;
    
    config->property_3D_value = NULL;
    config->binProperty1_3D_value = NULL;
    config->binProperty2_3D_value = NULL;
    config->binProperty3_3D_value = NULL;
    config->binProperty1_3D_mapLowLimit = NULL;
    config->binProperty1_3D_mapUpLimit = NULL;
    config->property_3D_value = NULL;
    config->binProperty1_3D_value = NULL;
    config->binProperty2_3D_value = NULL;
    config->binProperty3_3D_value = NULL;
    config->binsInLog1_3D_value = NULL;
    config->binsInLog2_3D_value = NULL;
    config->binsInLog3_3D_value = NULL;
    config->binsPerMag1_3D_value = NULL;
    config->binsPerMag2_3D_value = NULL;
    config->binsPerMag3_3D_value = NULL;    
    
    config->property_2D_value = NULL;
    config->binProperty1_2D_mapLowLimit = NULL;
    config->binProperty1_2D_mapUpLimit = NULL;
    config->binProperty1_2D_value = NULL;
    config->binProperty2_2D_value = NULL;
    config->binsInLog1_2D_value = NULL;
    config->binsInLog2_2D_value = NULL;
    config->binsPerMag1_2D_value = NULL;
    config->binsPerMag2_2D_value = NULL;

    config->property_1D_value = NULL;
    config->binProperty1_1D_mapLowLimit = NULL;
    config->binProperty1_1D_mapUpLimit = NULL;
    config->binProperty_1D_value = NULL;
    config->binsInLog_1D_value = NULL;
    config->binsPerMag_1D_value = NULL;    
    
    config->property_2D_median = NULL;
    config->binProperty1_2D_median_mapLowLimit = NULL;
    config->binProperty1_2D_median_mapUpLimit = NULL;
    config->binProperty1_2D_median = NULL;
    config->binProperty2_2D_median = NULL;
    config->binsInLog1_2D_median = NULL;
    config->binsInLog2_2D_median = NULL;
    config->binsPerMag1_2D_median = NULL;
    config->binsPerMag2_2D_median = NULL;
    
    config->binProperty1_2D = NULL;
    config->binProperty2_2D = NULL;
    config->binsInLog1_2D = NULL;
    config->binsInLog2_2D = NULL;
    config->binsPerMag1_2D = NULL;
    config->binsPerMag2_2D = NULL;
    
    config->binProperty_1D = NULL;
    config->binsInLog_1D = NULL;
    config->binsPerMag_1D = NULL;
    config->cumulative = NULL;
    
    config->outputLists = 0;
    
    config->property_1D_evolution = NULL;
    config->binProperty_1D_evolution = NULL;
    config->binsInLog_1D_evolution = NULL;
    config->binsPerMag_1D_evolution = NULL;
    config->binsMinValue_1D_evolution = NULL;
    config->binsMaxValue_1D_evolution = NULL;

    char *sps_model = NULL;
    
    //reading mandatory stuff
    getFromIni(&(config->inputType), parse_ini_get_int32,
               ini, "type", "Input");
    getFromIni(&(config->numFiles), parse_ini_get_int32,
               ini, "numFiles", "Input");
    getFromIni(&(config->inputFile), parse_ini_get_string,
               ini, "inputFile", "Input");
    getFromIni(&(config->boxsize), parse_ini_get_double,
               ini, "boxsize", "Input");
    getFromIni(&(config->gridsize), parse_ini_get_int32,
               ini, "gridsize", "Input");
    getFromIni(&(config->analytic), parse_ini_get_int32,
               ini, "analytic", "Input");
    getFromIni(&(config->mergertreeBinwidth), parse_ini_get_double,
               ini, "mergertreeBinwidthInLog", "Input");
    getFromIni(&(config->mergertreeEndRedshift), parse_ini_get_double,
               ini, "mergertreeEndRedshift", "Input");
    getFromIni(&(config->hmfFilename), parse_ini_get_string,
               ini, "hmfFiles", "Input");
    
    getFromIni(&(config->omega_m), parse_ini_get_double,
               ini, "omega_m", "Cosmology");
    getFromIni(&(config->omega_b), parse_ini_get_double,
               ini, "omega_b", "Cosmology");
    getFromIni(&(config->omega_l), parse_ini_get_double,
               ini, "omega_l", "Cosmology");
    getFromIni(&(config->hubble_h), parse_ini_get_double,
               ini, "hubble_h", "Cosmology");
     
    getFromIni(&(config->delayedSNfeedback), parse_ini_get_int32,
               ini, "doDelayedSNfeedback", "Simulation");
    getFromIni(&(config->fw), parse_ini_get_double,
               ini, "SNenergyFractionIntoWinds", "Simulation");
    getFromIni(&sps_model, parse_ini_get_string,
            ini, "stellarPopulationSynthesisModel", "Simulation");
    config->sps_model = 0;
    if(strcmp(sps_model, "S99") == 0)
    {
      config->sps_model = 1;
    }
    else if(strcmp(sps_model, "S99cont") == 0)
    {
      config->sps_model = 11;
    }
    else if(strcmp(sps_model, "BPASS") == 0)
    {
      config->sps_model = 2;
    }
    else if(strcmp(sps_model, "BPASScont") == 0)
    {
      config->sps_model = 12;
    }
    getFromIni(&(config->fesc), parse_ini_get_double,
               ini, "fesc", "Simulation");
    getFromIni(&(config->MvirThreshold), parse_ini_get_double,
               ini, "MvirThreshold", "Simulation");
    
    getFromIni(&(config->ionFilename), parse_ini_get_string,
               ini, "ionFilename", "Grid");
    getFromIni(&(config->densFilename), parse_ini_get_string,
               ini, "densFilename", "Grid");
    getFromIni(&(config->velxFilename), parse_ini_get_string,
               ini, "velxFilename", "Grid");
    getFromIni(&(config->velyFilename), parse_ini_get_string,
               ini, "velyFilename", "Grid");
    getFromIni(&(config->velzFilename), parse_ini_get_string,
               ini, "velzFilename", "Grid");
    getFromIni(&(config->ionInputInDoublePrecision), parse_ini_get_int32,
               ini, "ionInputInDoublePrecision", "Grid");
    getFromIni(&(config->densInputInDoublePrecision), parse_ini_get_int32,
               ini, "densInputInDoublePrecision", "Grid");
    getFromIni(&(config->velInputInDoublePrecision), parse_ini_get_int32,
               ini, "velInputInDoublePrecision", "Grid");
    getFromIni(&(config->memoryIntensive), parse_ini_get_int32,
               ini, "memoryIntensive", "Grid");
    getFromIni(&(config->smoothingScale), parse_ini_get_double,
               ini, "smoothingScale", "Grid");

    getFromIni(&(config->numOutputRedshifts), parse_ini_get_int32,
               ini, "numOutputRedshifts", "Analysis");
    if(config->numOutputRedshifts > 0)
    {
      getListFromIni(&(config->outputRedshifts), parse_ini_get_doublelist,
              ini, "outputRedshifts", "Analysis", config->numOutputRedshifts);
    }
    
    getFromIni(&(config->num_2D_history), parse_ini_get_int32,
               ini, "num_2D_history", "2dHistogram");
    if(config->num_2D_history > 0)
    {
      getListFromIni(&(config->property_2D_history), parse_ini_get_stringlist,
                ini, "property_2D_history", "2dHistogram", config->num_2D_history);
      getListFromIni(&(config->binProperty1_2D_history), parse_ini_get_stringlist,
                ini, "binProperty1_2D_history", "2dHistogram", config->num_2D_history);
      getListFromIni(&(config->binProperty2_2D_history), parse_ini_get_stringlist,
                ini, "binProperty2_2D_history", "2dHistogram", config->num_2D_history);
      getListFromIni(&(config->binsInLog1_2D_history), parse_ini_get_int32list,
                ini, "binsInLog1_2D_history", "2dHistogram", config->num_2D_history);
      getListFromIni(&(config->binsInLog2_2D_history), parse_ini_get_int32list,
                ini, "binsInLog2_2D_history", "2dHistogram", config->num_2D_history);
      getListFromIni(&(config->binsPerMag1_2D_history), parse_ini_get_int32list,
                ini, "binsPerMag1_2D_history", "2dHistogram", config->num_2D_history);
      getListFromIni(&(config->binsPerMag2_2D_history), parse_ini_get_int32list,
                ini, "binsPerMag2_2D_history", "2dHistogram", config->num_2D_history);
    }
    
    getFromIni(&(config->num_1D_history), parse_ini_get_int32,
               ini, "num_1D_history", "1dHistogram");
    if(config->num_1D_history > 0)
    {
      getListFromIni(&(config->property_1D_history), parse_ini_get_stringlist,
                ini, "property_1D_history", "1dHistogram", config->num_1D_history);
      getListFromIni(&(config->binProperty_1D_history), parse_ini_get_stringlist,
                ini, "binProperty_1D_history", "1dHistogram", config->num_1D_history);
      getListFromIni(&(config->binsInLog_1D_history), parse_ini_get_int32list,
                ini, "binsInLog_1D_history", "1dHistogram", config->num_1D_history);
      getListFromIni(&(config->binsPerMag_1D_history), parse_ini_get_int32list,
                ini, "binsPerMag_1D_history", "1dHistogram", config->num_1D_history);
    }
    
    getFromIni(&(config->num_2D_history_median), parse_ini_get_int32,
               ini, "num_2D_history_median", "2dHistogramHistoryMedian");
    if(config->num_2D_history_median > 0)
    {
      getListFromIni(&(config->property_2D_history_median), parse_ini_get_stringlist,
                ini, "property_2D_history_median", "2dHistogramHistoryMedian", config->num_2D_history_median);
      getListFromIni(&(config->binProperty1_2D_history_median), parse_ini_get_stringlist,
                ini, "binProperty1_2D_history_median", "2dHistogramHistoryMedian", config->num_2D_history_median);
      getListFromIni(&(config->binProperty2_2D_history_median), parse_ini_get_stringlist,
                ini, "binProperty2_2D_history_median", "2dHistogramHistoryMedian", config->num_2D_history_median);
      getListFromIni(&(config->binsInLog1_2D_history_median), parse_ini_get_int32list,
                ini, "binsInLog1_2D_history_median", "2dHistogramHistoryMedian", config->num_2D_history_median);
      getListFromIni(&(config->binsInLog2_2D_history_median), parse_ini_get_int32list,
                ini, "binsInLog2_2D_history_median", "2dHistogramHistoryMedian", config->num_2D_history_median);
      getListFromIni(&(config->binsPerMag1_2D_history_median), parse_ini_get_int32list,
                ini, "binsPerMag1_2D_history_median", "2dHistogramHistoryMedian", config->num_2D_history_median);
      getListFromIni(&(config->binsPerMag2_2D_history_median), parse_ini_get_int32list,
                ini, "binsPerMag2_2D_history_median", "2dHistogramHistoryMedian", config->num_2D_history_median);
    }
    
    getFromIni(&(config->num_3D_value), parse_ini_get_int32,
               ini, "num_3D_value", "3dHistogramValue");
    if(config->num_3D_value > 0)
    {
      getListFromIni(&(config->property_3D_value), parse_ini_get_stringlist,
                ini, "property_3D_value", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binProperty1_3D_mapLowLimit), parse_ini_get_doublelist,
                ini, "binProperty1_3D_mapLowLimit", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binProperty1_3D_mapUpLimit), parse_ini_get_doublelist,
                ini, "binProperty1_3D_mapUpLimit", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binProperty1_3D_value), parse_ini_get_stringlist,
                ini, "binProperty1_3D_value", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binProperty2_3D_value), parse_ini_get_stringlist,
                ini, "binProperty2_3D_value", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binProperty3_3D_value), parse_ini_get_stringlist,
                ini, "binProperty3_3D_value", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binsInLog1_3D_value), parse_ini_get_int32list,
                ini, "binsInLog1_3D_value", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binsInLog2_3D_value), parse_ini_get_int32list,
                ini, "binsInLog2_3D_value", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binsInLog3_3D_value), parse_ini_get_int32list,
                ini, "binsInLog3_3D_value", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binsPerMag1_3D_value), parse_ini_get_int32list,
                ini, "binsPerMag1_3D_value", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binsPerMag2_3D_value), parse_ini_get_int32list,
                ini, "binsPerMag2_3D_value", "3dHistogramValue", config->num_3D_value);
      getListFromIni(&(config->binsPerMag3_3D_value), parse_ini_get_int32list,
                ini, "binsPerMag3_3D_value", "3dHistogramValue", config->num_3D_value);
    }
    
    getFromIni(&(config->num_2D_value), parse_ini_get_int32,
               ini, "num_2D_value", "2dHistogramValue");
    if(config->num_2D_value > 0)
    {
      getListFromIni(&(config->property_2D_value), parse_ini_get_stringlist,
                ini, "property_2D_value", "2dHistogramValue", config->num_2D_value);
      getListFromIni(&(config->binProperty1_2D_mapLowLimit), parse_ini_get_doublelist,
                ini, "binProperty1_2D_mapLowLimit", "2dHistogramValue", config->num_2D_value);
      getListFromIni(&(config->binProperty1_2D_mapUpLimit), parse_ini_get_doublelist,
                ini, "binProperty1_2D_mapUpLimit", "2dHistogramValue", config->num_2D_value);
      getListFromIni(&(config->binProperty1_2D_value), parse_ini_get_stringlist,
                ini, "binProperty1_2D_value", "2dHistogramValue", config->num_2D_value);
      getListFromIni(&(config->binProperty2_2D_value), parse_ini_get_stringlist,
                ini, "binProperty2_2D_value", "2dHistogramValue", config->num_2D_value);
      getListFromIni(&(config->binsInLog1_2D_value), parse_ini_get_int32list,
                ini, "binsInLog1_2D_value", "2dHistogramValue", config->num_2D_value);
      getListFromIni(&(config->binsInLog2_2D_value), parse_ini_get_int32list,
                ini, "binsInLog2_2D_value", "2dHistogramValue", config->num_2D_value);
      getListFromIni(&(config->binsPerMag1_2D_value), parse_ini_get_int32list,
                ini, "binsPerMag1_2D_value", "2dHistogramValue", config->num_2D_value);
      getListFromIni(&(config->binsPerMag2_2D_value), parse_ini_get_int32list,
                ini, "binsPerMag2_2D_value", "2dHistogramValue", config->num_2D_value);
    }
    
    getFromIni(&(config->num_1D_value), parse_ini_get_int32,
               ini, "num_1D_value", "1dHistogramValue");
    if(config->num_1D_value > 0)
    {
      getListFromIni(&(config->property_1D_value), parse_ini_get_stringlist,
                ini, "property_1D_value", "1dHistogramValue", config->num_1D_value);
      getListFromIni(&(config->binProperty1_1D_mapLowLimit), parse_ini_get_doublelist,
                ini, "binProperty1_1D_mapLowLimit", "1dHistogramValue", config->num_1D_value);
      getListFromIni(&(config->binProperty1_1D_mapUpLimit), parse_ini_get_doublelist,
                ini, "binProperty1_1D_mapUpLimit", "1dHistogramValue", config->num_1D_value);
      getListFromIni(&(config->binProperty_1D_value), parse_ini_get_stringlist,
                ini, "binProperty_1D_value", "1dHistogramValue", config->num_1D_value);
      getListFromIni(&(config->binsInLog_1D_value), parse_ini_get_int32list,
                ini, "binsInLog_1D_value", "1dHistogramValue", config->num_1D_value);
      getListFromIni(&(config->binsPerMag_1D_value), parse_ini_get_int32list,
                ini, "binsPerMag_1D_value", "1dHistogramValue", config->num_1D_value);
    }    
    
    getFromIni(&(config->num_2D_median), parse_ini_get_int32,
               ini, "num_2D_median", "2dHistogramMedian");
    if(config->num_2D_median > 0)
    {
      getListFromIni(&(config->property_2D_median), parse_ini_get_stringlist,
                ini, "property_2D_median", "2dHistogramMedian", config->num_2D_median);
      getListFromIni(&(config->binProperty1_2D_median_mapLowLimit), parse_ini_get_doublelist,
                ini, "binProperty1_2D_median_mapLowLimit", "2dHistogramMedian", config->num_2D_median);
      getListFromIni(&(config->binProperty1_2D_median_mapUpLimit), parse_ini_get_doublelist,
                ini, "binProperty1_2D_median_mapUpLimit", "2dHistogramMedian", config->num_2D_median);
      getListFromIni(&(config->binProperty1_2D_median), parse_ini_get_stringlist,
                ini, "binProperty1_2D_median", "2dHistogramMedian", config->num_2D_median);
      getListFromIni(&(config->binProperty2_2D_median), parse_ini_get_stringlist,
                ini, "binProperty2_2D_median", "2dHistogramMedian", config->num_2D_median);
      getListFromIni(&(config->binsInLog1_2D_median), parse_ini_get_int32list,
                ini, "binsInLog1_2D_median", "2dHistogramMedian", config->num_2D_median);
      getListFromIni(&(config->binsInLog2_2D_median), parse_ini_get_int32list,
                ini, "binsInLog2_2D_median", "2dHistogramMedian", config->num_2D_median);
      getListFromIni(&(config->binsPerMag1_2D_median), parse_ini_get_int32list,
                ini, "binsPerMag1_2D_median", "2dHistogramMedian", config->num_2D_median);
      getListFromIni(&(config->binsPerMag2_2D_median), parse_ini_get_int32list,
                ini, "binsPerMag2_2D_median", "2dHistogramMedian", config->num_2D_median);
    }
    
    getFromIni(&(config->num_2D), parse_ini_get_int32,
               ini, "num_2D", "2dnumDensHistogram");
    if(config->num_2D > 0)
    {
      getListFromIni(&(config->binProperty1_2D), parse_ini_get_stringlist,
                ini, "binProperty1_2D", "2dnumDensHistogram", config->num_2D);
      getListFromIni(&(config->binProperty2_2D), parse_ini_get_stringlist,
                ini, "binProperty2_2D", "2dnumDensHistogram", config->num_2D);
      getListFromIni(&(config->binsInLog1_2D), parse_ini_get_int32list,
                ini, "binsInLog1_2D", "2dnumDensHistogram", config->num_2D);
      getListFromIni(&(config->binsInLog2_2D), parse_ini_get_int32list,
                ini, "binsInLog2_2D", "2dnumDensHistogram", config->num_2D);
      getListFromIni(&(config->binsPerMag1_2D), parse_ini_get_int32list,
                ini, "binsPerMag1_2D", "2dnumDensHistogram", config->num_2D);
      getListFromIni(&(config->binsPerMag2_2D), parse_ini_get_int32list,
                ini, "binsPerMag2_2D", "2dnumDensHistogram", config->num_2D);
    }
    
    getFromIni(&(config->num_1D), parse_ini_get_int32,
               ini, "num_1D", "1dnumDensHistogram");
    if(config->num_1D > 0)
    {
      getListFromIni(&(config->binProperty_1D), parse_ini_get_stringlist,
                ini, "binProperty_1D", "1dnumDensHistogram", config->num_1D);
      getListFromIni(&(config->binsInLog_1D), parse_ini_get_int32list,
                ini, "binsInLog_1D", "1dnumDensHistogram", config->num_1D);
      getListFromIni(&(config->binsPerMag_1D), parse_ini_get_int32list,
                ini, "binsPerMag_1D", "1dnumDensHistogram", config->num_1D);
      getListFromIni(&(config->cumulative), parse_ini_get_int32list,
                ini, "cumulative", "1dnumDensHistogram", config->num_1D);
    }
    
    getFromIni(&(config->trackEvolution), parse_ini_get_int32,
               ini, "trackEvolution", "AnalysisEvolution");
    
    getFromIni(&(config->num_1D_evolution), parse_ini_get_int32,
               ini, "num_1D_evolution", "1dHistogramEvolution");
    if(config->num_1D_evolution > 0)
    {
      getListFromIni(&(config->property_1D_evolution), parse_ini_get_stringlist,
                ini, "property_1D_evolution", "1dHistogramEvolution", config->num_1D_evolution);
      getListFromIni(&(config->binProperty_1D_evolution), parse_ini_get_stringlist,
                ini, "binProperty_1D_evolution", "1dHistogramEvolution", config->num_1D_evolution);
      getListFromIni(&(config->binsInLog_1D_evolution), parse_ini_get_int32list,
                ini, "binsInLog_1D_evolution", "1dHistogramEvolution", config->num_1D_evolution);
      getListFromIni(&(config->binsPerMag_1D_evolution), parse_ini_get_int32list,
                ini, "binsPerMag_1D_evolution", "1dHistogramEvolution", config->num_1D_evolution);
      getListFromIni(&(config->binsMinValue_1D_evolution), parse_ini_get_floatlist,
                ini, "binsMinValue_1D_evolution", "1dHistogramEvolution", config->num_1D_evolution);
      getListFromIni(&(config->binsMaxValue_1D_evolution), parse_ini_get_floatlist,
                ini, "binsMaxValue_1D_evolution", "1dHistogramEvolution", config->num_1D_evolution);
    }
    
    getFromIni(&(config->outputDir), parse_ini_get_string,
               ini, "outputDirectory", "Output");   
    getFromIni(&(config->outputLists), parse_ini_get_int32,
               ini, "writeTxtOutputLists", "Output");
    
    config->numSnaps = 0;
    config->scalefactors = NULL;
    config->redshifts = NULL;
    config->times = NULL;
    config->SNenergy = NULL;
    config->corrFactor_nion = NULL;
    
    config->size = 1;
    config->thisRank = 0;
    
    free(sps_model);
    
    return config;
}

extern void
dconfObj_del(dconfObj_t *config)
{
    assert(config != NULL);
    assert(*config != NULL);
    
    xfree((*config)->inputFile);
    xfree((*config)->hmfFilename);
    
    xfree((*config)->ionFilename);
    xfree((*config)->densFilename);
    xfree((*config)->velxFilename);
    xfree((*config)->velyFilename);
    xfree((*config)->velzFilename);

    xfree((*config)->outputRedshifts);
    
    for(int i=0; i<(*config)->num_2D_history; i++)
    {
      xfree((*config)->property_2D_history[i]);
      xfree((*config)->binProperty1_2D_history[i]);
      xfree((*config)->binProperty2_2D_history[i]);
    }
    xfree((*config)->property_2D_history);
    xfree((*config)->binProperty1_2D_history);
    xfree((*config)->binProperty2_2D_history);
    xfree((*config)->binsInLog1_2D_history);
    xfree((*config)->binsInLog2_2D_history);
    xfree((*config)->binsPerMag1_2D_history);
    xfree((*config)->binsPerMag2_2D_history);

    for(int i=0; i<(*config)->num_1D_history; i++)
    {
      xfree((*config)->property_1D_history[i]);
      xfree((*config)->binProperty_1D_history[i]);
    }
    xfree((*config)->property_1D_history);
    xfree((*config)->binProperty_1D_history);
    xfree((*config)->binsInLog_1D_history);
    xfree((*config)->binsPerMag_1D_history);
    
    
    for(int i=0; i<(*config)->num_2D_history_median; i++)
    {
      xfree((*config)->property_2D_history_median[i]);
      xfree((*config)->binProperty1_2D_history_median[i]);
      xfree((*config)->binProperty2_2D_history_median[i]);
    }
    xfree((*config)->property_2D_history_median);
    xfree((*config)->binProperty1_2D_history_median);
    xfree((*config)->binProperty2_2D_history_median);
    xfree((*config)->binsInLog1_2D_history_median);
    xfree((*config)->binsInLog2_2D_history_median);
    xfree((*config)->binsPerMag1_2D_history_median);
    xfree((*config)->binsPerMag2_2D_history_median);
    
    
    for(int i=0; i<(*config)->num_3D_value; i++)
    {
      xfree((*config)->property_3D_value[i]);
      xfree((*config)->binProperty1_3D_value[i]);
      xfree((*config)->binProperty2_3D_value[i]);
      xfree((*config)->binProperty3_3D_value[i]);
    }
    xfree((*config)->binProperty1_3D_mapLowLimit);
    xfree((*config)->binProperty1_3D_mapUpLimit);
    xfree((*config)->property_3D_value);
    xfree((*config)->binProperty1_3D_value);
    xfree((*config)->binProperty2_3D_value);
    xfree((*config)->binProperty3_3D_value);
    xfree((*config)->binsInLog1_3D_value);
    xfree((*config)->binsInLog2_3D_value);
    xfree((*config)->binsInLog3_3D_value);
    xfree((*config)->binsPerMag1_3D_value);
    xfree((*config)->binsPerMag2_3D_value);
    xfree((*config)->binsPerMag3_3D_value);
    
    for(int i=0; i<(*config)->num_2D_value; i++)
    {
      xfree((*config)->property_2D_value[i]);
      xfree((*config)->binProperty1_2D_value[i]);
      xfree((*config)->binProperty2_2D_value[i]);
    }
    xfree((*config)->binProperty1_2D_mapLowLimit);
    xfree((*config)->binProperty1_2D_mapUpLimit);
    xfree((*config)->property_2D_value);
    xfree((*config)->binProperty1_2D_value);
    xfree((*config)->binProperty2_2D_value);
    xfree((*config)->binsInLog1_2D_value);
    xfree((*config)->binsInLog2_2D_value);
    xfree((*config)->binsPerMag1_2D_value);
    xfree((*config)->binsPerMag2_2D_value);

    for(int i=0; i<(*config)->num_1D_value; i++)
    {
      xfree((*config)->property_1D_value[i]);
      xfree((*config)->binProperty_1D_value[i]);
    }
    xfree((*config)->binProperty1_1D_mapLowLimit);
    xfree((*config)->binProperty1_1D_mapUpLimit);
    xfree((*config)->property_1D_value);
    xfree((*config)->binProperty_1D_value);
    xfree((*config)->binsInLog_1D_value);
    xfree((*config)->binsPerMag_1D_value);
    
    for(int i=0; i<(*config)->num_2D_median; i++)
    {
      xfree((*config)->property_2D_median[i]);
      xfree((*config)->binProperty1_2D_median[i]);
      xfree((*config)->binProperty2_2D_median[i]);
    }
    xfree((*config)->binProperty1_2D_median_mapLowLimit);
    xfree((*config)->binProperty1_2D_median_mapUpLimit);
    xfree((*config)->property_2D_median);
    xfree((*config)->binProperty1_2D_median);
    xfree((*config)->binProperty2_2D_median);
    xfree((*config)->binsInLog1_2D_median);
    xfree((*config)->binsInLog2_2D_median);
    xfree((*config)->binsPerMag1_2D_median);
    xfree((*config)->binsPerMag2_2D_median);    
    
    for(int i=0; i<(*config)->num_2D; i++)
    {
      xfree((*config)->binProperty1_2D[i]);
      xfree((*config)->binProperty2_2D[i]);
    }
    xfree((*config)->binProperty1_2D);
    xfree((*config)->binProperty2_2D);
    xfree((*config)->binsInLog1_2D);
    xfree((*config)->binsInLog2_2D);
    xfree((*config)->binsPerMag1_2D);
    xfree((*config)->binsPerMag2_2D);
    
    for(int i=0; i<(*config)->num_1D; i++)
    {
      xfree((*config)->binProperty_1D[i]);
    }
    xfree((*config)->binProperty_1D);
    xfree((*config)->binsInLog_1D);
    xfree((*config)->binsPerMag_1D);
    xfree((*config)->cumulative);
    
    for(int i=0; i<(*config)->num_1D_evolution; i++)
    {
      xfree((*config)->property_1D_evolution[i]);
      xfree((*config)->binProperty_1D_evolution[i]);
    }
    xfree((*config)->property_1D_evolution);
    xfree((*config)->binProperty_1D_evolution);
    xfree((*config)->binsInLog_1D_evolution);
    xfree((*config)->binsPerMag_1D_evolution);
    xfree((*config)->binsMinValue_1D_evolution);
    xfree((*config)->binsMaxValue_1D_evolution);

    xfree((*config)->outputDir);
    
    xfree((*config)->scalefactors);
    xfree((*config)->redshifts);
    xfree((*config)->times);
    if((*config)->SNenergy != NULL) xfree((*config)->SNenergy);
    if((*config)->corrFactor_nion != NULL) xfree((*config)->corrFactor_nion);

    xfree(*config);
    *config = NULL;
}

extern dconfObj_t
readDconfObj(char *fileName) {
    dconfObj_t newConfObj;
    parse_ini_t ini;
    
    ini = parse_ini_open(fileName);
    if (ini == NULL) {
        fprintf(stderr, "FATAL:  Could not open %s for reading.\n",
                fileName);
        exit(EXIT_FAILURE);
    }
    
    newConfObj = dconfObj_new(ini);
    
    parse_ini_close(&ini);
    
    return newConfObj;
}

/*--- Implementations of local functions --------------------------------*/
