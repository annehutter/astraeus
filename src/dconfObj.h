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
    
    double         omega_m;
    double         omega_b;
    double         omega_l;
    double         hubble_h;
    
    double         FS;
    float          *timeSnaps;
    
    int            delayedSNfeedback;
    float          *SNenergy;
    double         FW;
    
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
