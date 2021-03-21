/*
 *  confObj.h
 *  uvff
 *
 *  Created by 
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef DELPHI_CONFOBJ_H
#define DELPHI_CONFOBJ_H

/*--- Includes ----------------------------------------------------------*/
#include "parse_ini.h"
#include <stdint.h>
#include <stdbool.h>


/*--- ADT handle --------------------------------------------------------*/
typedef struct dconfObj_struct *dconfObj_t;


/*--- Implemention of main structure ------------------------------------*/
struct dconfObj_struct {
    int         inputType;
    int         numFiles;
    char        *inputFile;
    double      boxsize;
    int         gridsize;
    int         analytic;
    double      mergertreeBinwidth;
    double      mergertreeEndRedshift;
    char        *hmfFilename;
    
    double      omega_m;
    double      omega_b;
    double      omega_l;
    double      hubble_h;
    
    int         delayedSNfeedback;
    double      fw;
    int         sps_model;
    double      fesc;
    double      MvirThreshold;
    
    char        *ionFilename;
    char        *densFilename;
    char        *velxFilename;
    char        *velyFilename;
    char        *velzFilename;
    int         ionInputInDoublePrecision;
    int         densInputInDoublePrecision;
    int         velInputInDoublePrecision;
    int         memoryIntensive;
    double      smoothingScale;
    
    int         numOutputRedshifts;
    double      *outputRedshifts;
    
    int         num_2D_history;
    char        **property_2D_history;
    char        **binProperty1_2D_history;
    char        **binProperty2_2D_history;
    int         *binsInLog1_2D_history;
    int         *binsInLog2_2D_history;
    int         *binsPerMag1_2D_history;
    int         *binsPerMag2_2D_history;

    int         num_1D_history;
    char        **property_1D_history;
    char        **binProperty_1D_history;
    int         *binsInLog_1D_history;
    int         *binsPerMag_1D_history;
    
    int         num_2D_history_median;
    char        **property_2D_history_median;
    char        **binProperty1_2D_history_median;
    char        **binProperty2_2D_history_median;
    int         *binsInLog1_2D_history_median;
    int         *binsInLog2_2D_history_median;
    int         *binsPerMag1_2D_history_median;
    int         *binsPerMag2_2D_history_median;
    
    int         num_3D_value;
    char        **property_3D_value;
    double      *binProperty1_3D_mapLowLimit;
    double      *binProperty1_3D_mapUpLimit;
    char        **binProperty1_3D_value;
    char        **binProperty2_3D_value;
    char        **binProperty3_3D_value;
    int         *binsInLog1_3D_value;
    int         *binsInLog2_3D_value;
    int         *binsInLog3_3D_value;
    int         *binsPerMag1_3D_value;
    int         *binsPerMag2_3D_value;
    int         *binsPerMag3_3D_value;
    
    int         num_2D_value;
    char        **property_2D_value;
    double      *binProperty1_2D_mapLowLimit;
    double      *binProperty1_2D_mapUpLimit;
    char        **binProperty1_2D_value;
    char        **binProperty2_2D_value;
    int         *binsInLog1_2D_value;
    int         *binsInLog2_2D_value;
    int         *binsPerMag1_2D_value;
    int         *binsPerMag2_2D_value;

    int         num_1D_value;
    double      *binProperty1_1D_mapLowLimit;
    double      *binProperty1_1D_mapUpLimit;
    char        **property_1D_value;
    char        **binProperty_1D_value;
    int         *binsInLog_1D_value;
    int         *binsPerMag_1D_value;
    
    int         num_2D_median;
    char        **property_2D_median;
    double      *binProperty1_2D_median_mapLowLimit;
    double      *binProperty1_2D_median_mapUpLimit;
    char        **binProperty1_2D_median;
    char        **binProperty2_2D_median;
    int         *binsInLog1_2D_median;
    int         *binsInLog2_2D_median;
    int         *binsPerMag1_2D_median;
    int         *binsPerMag2_2D_median;
    
    int         num_2D;
    char        **binProperty1_2D;
    char        **binProperty2_2D;
    int         *binsInLog1_2D;
    int         *binsInLog2_2D;
    int         *binsPerMag1_2D;
    int         *binsPerMag2_2D;
    
    int         num_1D;
    char        **binProperty_1D;
    int         *binsInLog_1D;
    int         *binsPerMag_1D;
    int         *cumulative;
        
    int         trackEvolution;
    
    int         num_1D_evolution;
    char        **property_1D_evolution;
    char        **binProperty_1D_evolution;
    int         *binsInLog_1D_evolution;
    int         *binsPerMag_1D_evolution;
    float       *binsMinValue_1D_evolution;
    float       *binsMaxValue_1D_evolution;
    
    char        *outputDir;
    int         outputLists;

    /* internal variables */
    int         numSnaps;
    float       *scalefactors;
    float       *redshifts;
    float       *times;
    float       *SNenergy;
    double      *corrFactor_nion;
    
    int         thisRank;
    int         size;
};


/*--- Prototypes of exported functions ----------------------------------*/
extern dconfObj_t
readDconfObj(char *fileName);

extern dconfObj_t
dconfObj_new(parse_ini_t ini);

extern void
dconfObj_del(dconfObj_t *config);


#endif
