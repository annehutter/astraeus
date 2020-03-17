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
    
    config->outputRedshifts = NULL;
    
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
    
    config->property_1D_evolution = NULL;
    config->binProperty_1D_evolution = NULL;
    config->binsInLog_1D_evolution = NULL;
    config->binsPerMag_1D_evolution = NULL;
    config->binsMinValue_1D_evolution = NULL;
    config->binsMaxValue_1D_evolution = NULL;

    //reading mandatory stuff
    getFromIni(&(config->numFiles), parse_ini_get_int32,
               ini, "numFiles", "Input");
    getFromIni(&(config->inputFile), parse_ini_get_string,
               ini, "inputFile", "Input");
    getFromIni(&(config->boxsize), parse_ini_get_double,
               ini, "boxsize", "Input");
    getFromIni(&(config->gridsize), parse_ini_get_int32,
               ini, "gridsize", "Input");
    
    getFromIni(&(config->omega_m), parse_ini_get_double,
               ini, "omega_m", "Cosmology");
    getFromIni(&(config->omega_b), parse_ini_get_double,
               ini, "omega_b", "Cosmology");
    getFromIni(&(config->omega_l), parse_ini_get_double,
               ini, "omega_l", "Cosmology");
    getFromIni(&(config->hubble_h), parse_ini_get_double,
               ini, "hubble_h", "Cosmology");
     
    getFromIni(&(config->ionFilename), parse_ini_get_string,
               ini, "ionFilename", "Grid");
    getFromIni(&(config->densFilename), parse_ini_get_string,
               ini, "densFilename", "Grid");
    getFromIni(&(config->ionInputInDoublePrecision), parse_ini_get_int32,
               ini, "ionInputInDoublePrecision", "Grid");
    getFromIni(&(config->densInputInDoublePrecision), parse_ini_get_int32,
               ini, "densInputInDoublePrecision", "Grid");
    getFromIni(&(config->memoryIntensive), parse_ini_get_int32,
               ini, "memoryIntensive", "Grid");
    
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
    
    config->numSnaps = 0;
    config->scalefactors = NULL;
    config->redshifts = NULL;
    config->times = NULL;
    
    config->size = 1;
    config->thisRank = 0;
    
    return config;
}

extern void
dconfObj_del(dconfObj_t *config)
{
    assert(config != NULL);
    assert(*config != NULL);
    
    xfree((*config)->inputFile);
    
    xfree((*config)->ionFilename);
    xfree((*config)->densFilename);

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
