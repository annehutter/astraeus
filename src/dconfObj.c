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


/*--- Implementations of exported functios ------------------------------*/
extern dconfObj_t
dconfObj_new(parse_ini_t ini)
{
    dconfObj_t config;
    assert(ini != NULL);
    
    config = xmalloc(sizeof(struct dconfObj_struct));
    
    char *reion_model = NULL;
    char *sps_model = NULL;
    char *fesc_model = NULL;
    char *radfeedback_model = NULL;
    
    //reading mandatory stuff
    getFromIni(&(config->numFiles), parse_ini_get_int32,
               ini, "numFiles", "Input");
    getFromIni(&(config->fileName), parse_ini_get_string,
               ini, "fileName", "Input");
    
    getFromIni(&(config->redshiftFile), parse_ini_get_string,
               ini, "redshiftFile", "Simulation");
    getFromIni(&(config->startSnap), parse_ini_get_int32,
               ini, "startSnapshot", "Simulation");
    getFromIni(&(config->endSnap), parse_ini_get_int32,
               ini, "endSnapshot", "Simulation");
    getFromIni(&(config->deltaSnap), parse_ini_get_int32,
               ini, "deltaSnapshot", "Simulation");
        
    getFromIni(&(config->gridsize), parse_ini_get_int32,
               ini, "gridsize", "Simulation");
    getFromIni(&(config->boxsize), parse_ini_get_double,
               ini, "boxsize", "Simulation");
    config->inv_boxsize = 1./config->boxsize;
    
    getFromIni(&(config->memoryIntensive), parse_ini_get_int32,
               ini, "fastButMemoryIntensive", "Simulation");
    
    getFromIni(&(config->omega_m), parse_ini_get_double,
               ini, "OM0", "Cosmology");
    getFromIni(&(config->omega_b), parse_ini_get_double,
               ini, "OB0", "Cosmology");
    getFromIni(&(config->omega_l), parse_ini_get_double,
               ini, "OL0", "Cosmology");
    getFromIni(&(config->hubble_h), parse_ini_get_double,
               ini, "HUBBLE_CONSTANT", "Cosmology");
    
    config->dt_model = 0;
    config->dt_rescaleFactor = 1.;
    config->dt_deltaTimeInMyr = 0.;
    checkFromIni(&(config->dt_model), parse_ini_get_int32,
               ini, "timestepModel", "StarFormation");
    if(config->dt_model == 1)
    {
      getFromIni(&(config->dt_rescaleFactor), parse_ini_get_float,
               ini, "timestepModel1_factorRescale", "StarFormation");
    }
    if(config->dt_model > 1)
    {
      getFromIni(&(config->dt_deltaTimeInMyr), parse_ini_get_float,
               ini, "timestepModel2_deltaTimeInMyr", "StarFormation");
    }
    
    getFromIni(&(config->FS), parse_ini_get_double,
               ini, "starFormationEfficiency", "StarFormation");
    config->timeSnaps = NULL;
    
    getFromIni(&(config->delayedSNfeedback), parse_ini_get_int32,
               ini, "doDelayedSNfeedback", "SNfeedback");
    config->SNenergy = NULL;
    getFromIni(&(config->FW), parse_ini_get_double,
               ini, "SNenergyFractionIntoWinds", "SNfeedback");
    
    getFromIni(&(config->radfeedback), parse_ini_get_int32,
               ini, "doRadfeedback", "RadiativeFeedback");
    getFromIni(&radfeedback_model, parse_ini_get_string,
               ini, "radfeedbackModel", "RadiativeFeedback");
    config->radfeedback_model = -1;
    if(config->radfeedback == 1)
    {
      if(strcmp(radfeedback_model, "SOBACCHI") == 0)
      {
        config->radfeedback_model = 1;
      }
      else if(strcmp(radfeedback_model, "GNEDIN") == 0)
      {
        config->radfeedback_model = 2;
      }
      else if(strcmp(radfeedback_model, "MJEANS") == 0)
      {
        config->radfeedback_model = 3;
      }
      else if(strcmp(radfeedback_model, "MIN") == 0)
      {
        config->radfeedback_model = 4;
      }
      else if(strcmp(radfeedback_model, "TEMPEVOL") == 0)
      {
        config->radfeedback_model = 5;
      }
      else
      {
        fprintf(stderr, "This reionization model is not supported.\n");
        exit(EXIT_FAILURE);
      }
    }
    getFromIni(&(config->XHII_threshold), parse_ini_get_double,
               ini, "ionThreshold", "RadiativeFeedback");
    getFromIni(&(config->temp), parse_ini_get_float,
               ini, "tempIonGas", "RadiativeFeedback");
    getFromIni(&(config->mu), parse_ini_get_float,
               ini, "muGas", "RadiativeFeedback");
     
    getFromIni(&(config->reion), parse_ini_get_int32,
               ini, "doReionization", "Reionization");
    getFromIni(&config->cifogIniFile, parse_ini_get_string,
               ini, "cifogIniFile", "Reionization");
    
    getFromIni(&reion_model, parse_ini_get_string,
               ini, "reionizationModel", "Reionization");
    config->reion_model = -1;
    if(config->reion == 1)
    {
      if(strcmp(reion_model, "LOCAL") == 0)
      {
        config->reion_model = 1;
      }
      else
      {
        fprintf(stderr, "This reionization model is not supported.\n");
        exit(EXIT_FAILURE);
      }
    }
    config->corrFactor_nion = NULL;
    
    getFromIni(&sps_model, parse_ini_get_string,
            ini, "stellarPopulationSynthesisModel", "Reionization");
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
    
    getFromIni(&fesc_model, parse_ini_get_string,
               ini, "fescModel", "Reionization");   
    config->fesc_model = -1;
    if(config->reion == 1)
    {
      if(strcmp(fesc_model, "MHDEC") == 0)
      {
        config->fesc_model = 1;
      }
      else if(strcmp(fesc_model, "MHINC") == 0)
      {
        config->fesc_model = 2;
      }
      else if(strcmp(fesc_model, "CONST") == 0)
      {
        config->fesc_model = 0;
      }
      else if(strcmp(fesc_model, "SN") == 0)
      {
        config->fesc_model = 3;
      }
      else
      {
        fprintf(stderr, "This fesc model is not supported.\n");
        exit(EXIT_FAILURE);
      }
    }
    config->fesc = 1.;
    config->MHlow = 1.e8;
    config->MHhigh = 1.e12;
    config->fescLow = 0.99;
    config->fescHigh = 0.1;
    if(config->fesc_model == 1 || config->fesc_model == 2)
    {
      getFromIni(&(config->MHlow), parse_ini_get_double,
               ini, "MHlow", "fescMH");
      getFromIni(&(config->MHhigh), parse_ini_get_double,
               ini, "MHhigh", "fescMH");
      getFromIni(&(config->fescLow), parse_ini_get_double,
               ini, "fescLow", "fescMH");
      getFromIni(&(config->fescHigh), parse_ini_get_double,
               ini, "fescHigh", "fescMH");
    }
    else if(config->fesc_model == 0)
    {
      getFromIni(&(config->fesc), parse_ini_get_double,
               ini, "fesc", "fescConst");
    }
    else if(config->fesc_model == 3)
    {
      getFromIni(&(config->fesc), parse_ini_get_double,
               ini, "fesc", "fescSN");
    }
    
    config->outputType = 1;
    getFromIni(&(config->outputType), parse_ini_get_int32,
               ini, "type", "Output");
    getFromIni(&(config->numSnapsToWrite), parse_ini_get_int32,
               ini, "numSnapsToWrite", "Output");
    getFromIni(&(config->horizontalOutput), parse_ini_get_int32,
               ini, "horizontalOutput", "Output");
               
    if((config->numSnapsToWrite > 0)&&(config->horizontalOutput==1))
    {
      getListFromIni(&(config->snapList), parse_ini_get_int32list,
                ini, "snapList", "Output", config->numSnapsToWrite);
    }
    else
    {
      config->snapList = xmalloc(sizeof(int32_t *));
      config->snapList[0] = -1;
    }
    
    getFromIni(&(config->verticalOutput), parse_ini_get_int32,
               ini, "verticalOutput", "Output");
    getFromIni(&(config->percentageOfTreesToWrite), parse_ini_get_int32,
               ini, "percentageOfTreesToWrite", "Output");
               
    getFromIni(&(config->outFileName), parse_ini_get_string,
                ini, "outputFile", "Output");
    if(directory_exist(config->outFileName) != 1)
    {
      printf("Directory for output does not exist. This will result in a segmentation fault at the end of the program.\nAborting now!\n");
      exit(EXIT_FAILURE);
    }
               
    xfree(reion_model);
    xfree(sps_model);
    xfree(fesc_model);
    xfree(radfeedback_model);

    return config;
}

extern void
dconfObj_del(dconfObj_t *config)
{
    assert(config != NULL);
    assert(*config != NULL);
    
    //General
    xfree((*config)->fileName);
    xfree((*config)->snapList);
    xfree((*config)->redshiftFile);
    if((*config)->timeSnaps != NULL) xfree((*config)->timeSnaps);
    if((*config)->SNenergy != NULL) xfree((*config)->SNenergy);
    if((*config)->corrFactor_nion != NULL) xfree((*config)->corrFactor_nion);
    xfree((*config)->cifogIniFile);
    xfree((*config)->outFileName);

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
