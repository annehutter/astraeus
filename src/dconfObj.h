/*
 *  confObj.h
 *  uvff
 *
 *  Created by
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DELPHI_CONFOBJ_H
#define DELPGI_CONFOBJ_H

/*--- Includes ----------------------------------------------------------*/
#include "parse_ini.h"
#include <stdint.h>
#include <stdbool.h>


/*--- ADT handle --------------------------------------------------------*/
typedef struct dconfObj_struct *dconfObj_t;


/*--- Implemention of main structure ------------------------------------*/
struct dconfObj_struct {
    //Input SAGE galaxies
    int            numFiles;
    char           *fileName;

    char           *redshiftFile;
    int            startSnap;
    int            endSnap;
    int            deltaSnap;
    int            gridsize;
    double         boxsize;
    double         inv_boxsize;
    int            memoryIntensive;

    double         omega_m;
    double         omega_b;
    double         omega_l;
    double         hubble_h;
    
    int            dt_model;
    float          dt_rescaleFactor;
    float          dt_deltaTimeInMyr;

    double         FS;
    float          *timeSnaps;

    int            delayedSNfeedback;
    float          *SNenergy;
    double         FW;
    int            burstySF;

#if defined WITHMETALS
    /* For metals and dust */
    int            metals;
    char           *metal_table_file;
    double         metal_ejectLoadingFactor;

    int            metalTables_numMassBins;
    int            metalTables_numMetallicityBins;
    double         *metalTables_mass;
    double         *metalTables_metallicity;
    double         *metalTables_imf;
    double         *metalTables_lifetime;
    double         *metalTables_metalYields;
    double         *metalTables_metalYields_SNIa;
    double         *metalTables_metalYieldsO;
    double         *metalTables_metalYieldsO_SNIa;
    double         *metalTables_metalYieldsFe;
    double         *metalTables_metalYieldsFe_SNIa;
    double         *metalTables_gasYields;
    double         *metalTables_dustYields;
    /* invdeltat is useful to compute the SFR */
    double         *invdeltat;
    
    double         dust_yield_SNII;
    double         fracColdGas;
    double         dust_destrEfficiency;
    double         dust_timescaleGrainGrowth;
#endif
    
    int            radfeedback;
    int            radfeedback_model;
    double         XHII_threshold;
    float          temp;
    float          mu;

    int            reion;
    char           *cifogIniFile;
    int            reion_model;
    int            sps_model;
    double         *corrFactor_nion;
    int            fesc_model;
    double         fesc;
    double         MHlow;
    double         MHhigh;
    double         fescLow;
    double         fescHigh;

    int            outputType;
    int            horizontalOutput;
    int            numSnapsToWrite;
    int            *snapList;
    int            verticalOutput;
    int            percentageOfTreesToWrite;
    char           *outFileName;
};


/*--- Prototypes of exported functions ----------------------------------*/
extern dconfObj_t
readDconfObj(char *fileName);

extern dconfObj_t
dconfObj_new(parse_ini_t ini);

extern void
dconfObj_del(dconfObj_t *config);


#endif
