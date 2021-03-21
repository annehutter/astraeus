#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <stdint.h>

#ifdef __MPI
#include <mpi.h>
#endif

#include "phys_const.h"
#include "utils.h"
#include "dconfObj.h"

/* NEEDED: times from chosen snapshot to past */
/* UV luminosity for those times */

float convert_UVperFreq_to_UVmag(float UVperFreq)
{
    float UVmag = -2.5 * log10(UVperFreq / (4.*PI*100.*pc_cm*pc_cm)) - 48.6;
    
    return UVmag;
}

float get_UVperFreq_after_time(float timeInYr)
{
    float logUVperAngstroem = 33.44 - 1.33*0.30;
    if(timeInYr >= 4.e6)
        logUVperAngstroem = -1.33 * log10(timeInYr * 5.e-7) + 33.44;
    
    float UVperFreq = pow(10., logUVperAngstroem) * lambda_UV_cm*lambda_UV_cm*1.e8 / clight_cm;
    
    return UVperFreq;
}

float get_UVcorrectionFactor(float timeOffset, float timeBegin, float timeEnd, float timeBreak)
{
    float deltaTime = timeEnd - timeBegin;
    
    float UVcorrectionFactor = 1.;
    
    if(deltaTime>4.e6)
    {
        UVcorrectionFactor = timeBreak*(1. - 3.*(pow(deltaTime*2.5e-7, -0.33)-1.))/deltaTime;
    }
 
    if(timeBegin-timeOffset>=4.e6)
    {
        UVcorrectionFactor = - 3.*(timeBegin-timeOffset)*(pow((timeEnd-timeOffset)/(timeBegin-timeOffset), -0.33) - 1.)/deltaTime;
    }
    
    return UVcorrectionFactor;
}

float calc_UV_lum_SMH(int numSnaps, double *starmassHistoryAll, int thisGal, int currSnap, float *corrFactor)
{
    double *starmassHistory = &(starmassHistoryAll[(long int)thisGal*numSnaps]);
    
    /* sum up contribution from previous times */
    float totUVperFreq = 0.;

    for(long int i=0; i<=currSnap; i++)
    {
        totUVperFreq += corrFactor[i] * starmassHistory[i];
    }
        
    return totUVperFreq;
}

float *calc_UV_lum_at_snap(dconfObj_t simParam, int numGal, int numSnaps, double *starmassHistoryAll, int currSnap)
{
    float *snaps_ages = simParam->times;
    
    /* get time of current snapshot */
    float currTime = snaps_ages[currSnap];
    
    /* get differences to previous times */
    int numPrevTimes = currSnap + 1;
    float *prevTimes = allocate_array_float(numPrevTimes , "previousTimes");
    float *corrFactor = allocate_array_float(numPrevTimes, "corrFactor");
    
    for(int prevSnap=0; prevSnap<=currSnap; prevSnap++)
        prevTimes[prevSnap] = (currTime - snaps_ages[prevSnap]) / yr_s;
        
    /* correct by factor that stellar mass forms continuously in the time bin */
    corrFactor[0] = get_UVperFreq_after_time(prevTimes[0]);
    for(int prevSnap=1; prevSnap<=currSnap; prevSnap++) 
    {
        corrFactor[prevSnap] = get_UVcorrectionFactor(prevTimes[currSnap], prevTimes[prevSnap], prevTimes[prevSnap-1], 4.e6) * get_UVperFreq_after_time(prevTimes[prevSnap]);
    }
 
    float *Lc_int = allocate_array_float(numGal, "Lc_int");
        
    for(int gal=0; gal<numGal; gal++)
    {
        Lc_int[gal] = calc_UV_lum_SMH(numSnaps, starmassHistoryAll, gal, currSnap, corrFactor);
    }
    
    free(corrFactor);
    free(prevTimes);
    
    return Lc_int;
}

double *calc_UV_mag(dconfObj_t simParam, int numGal, int numSnaps, double *starmassHistoryAll, int currSnap)
{    
    float *Lc = calc_UV_lum_at_snap(simParam, numGal, numSnaps, starmassHistoryAll, currSnap);
    double *UVmag = allocate_array_double(numGal, "UVmag");
    
    /* Convert total UVperFreq into magnitudes */
    for(int gal=0; gal<numGal; gal++)
    {
        if(Lc[gal] > 0.)
            UVmag[gal] = convert_UVperFreq_to_UVmag(Lc[gal]);
        else
            UVmag[gal] = 0.;
    }
    
    free(Lc);
    
    return UVmag;
}