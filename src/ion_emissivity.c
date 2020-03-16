#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include "utils.h"
#include "dconfObj.h"
#include "gal_gtree.h"
#include "fesc.h"
#include "ion_emissivity.h"

double get_nion_for_model(gal_t *thisGal, dconfObj_t simParam)
{
  double nion_sps = get_nion_sps(thisGal, simParam);
  double fesc = get_fesc(thisGal, simParam);
    
  return fesc * nion_sps;
}

double get_nion_sps(gal_t *thisGal, dconfObj_t simParam)
{
  float *times = simParam->timeSnaps;
  int sps_model = simParam->sps_model;
  double *corrFactor_nion = simParam->corrFactor_nion;
  double nion_sps = 0.;
  
  if(sps_model == 1)
    nion_sps = get_nion_S99(thisGal, times)/simParam->hubble_h;
  else if(sps_model == 2)
    nion_sps = get_nion_BPASS(thisGal, times)/simParam->hubble_h;
  else if(sps_model == 12 || sps_model == 11)
    nion_sps = get_nion_cont(thisGal, corrFactor_nion)/simParam->hubble_h;
  else
    nion_sps = thisGal->Mvir * 1.5e43;
  
  return nion_sps;
}

double get_nion_cont(gal_t *thisGal, double *corrFactor_nion)
{ 
  int snapGal = thisGal->snapnumber;
  double nion = 0.;
  
  for(int snap=1; snap<=snapGal; snap++)
  {
    nion += thisGal->stellarmasshistory[snap] * corrFactor_nion[snapGal*(snapGal-1)/2 + snap];
  }
  
  return nion;
}

double get_nion_S99(gal_t *thisGal, float *times)
{ 
  int snapGal = thisGal->snapnumber;
  float age = 0.;
  double nion = 0.;
  
  for(int snap=1; snap<=snapGal; snap++)
  {
    age = times[snapGal] - times[snap];
    nion += thisGal->stellarmasshistory[snap] * nion_S99(age);
  }
  
  return nion;
}

double nion_S99(float age)
{
  double MyrInSec = 3.1536e13;
  double MyrInSec_inv = 3.1709792e-14;

  if(age < 3.16*MyrInSec)
    age = 3.16*MyrInSec;
  
  double nion = 2.18e47 * pow(0.5 * age * MyrInSec_inv, -3.92);

  return nion;
}

double get_nion_BPASS(gal_t *thisGal, float *times)
{
  int snapGal = thisGal->snapnumber;
  float age = 0.;
  double nion = 0.;
  
  for(int snap=1; snap<=snapGal; snap++)
  {
    age = times[snapGal] - times[snap];
    nion += thisGal->stellarmasshistory[snap] * nion_BPASS(age);
  }
  
  return nion;
}

double nion_BPASS(float age)
{  
  double MyrInSec = 3.1536e13;
  double MyrInSec_inv = 3.1709792e-14;
  
  if(age < 3.16*MyrInSec)
    age = 3.16*MyrInSec;
  
  double nion = 9.091836e46 * pow(0.5 * age * MyrInSec_inv, -2.28);

  return nion;
}


/*---------------------------------------------------------------------*/
/* LOOKUP TABLE                                                        */
/*---------------------------------------------------------------------*/

double *get_corrFactor_nion(dconfObj_t simParam)
{
  int endSnap = simParam->endSnap;
  float *times = simParam->timeSnaps;
  int sps_model = simParam->sps_model;
  double *corrFactor_nion = NULL;
  
  if(sps_model == 11)
  {
    corrFactor_nion = get_corrFactor_nion_S99(endSnap, times);
  }
  else if(sps_model == 12)
  {
    corrFactor_nion = get_corrFactor_nion_BPASS(endSnap, times);
  }
  else
  {
    fprintf(stderr, "Not supported stellar population synthesis model for number of ionizing photons. Please revise!\n");
    exit(EXIT_FAILURE);
  }
  
  return corrFactor_nion;
}

double *get_corrFactor_nion_S99(int endSnap, float *times)
{
  float MyrInSec = 3.1536e13;
  float t0 = 2. * MyrInSec;
  float tbreak = 3.16 * MyrInSec;
  float exponent = -3.92;
  double nion = 2.183735e47;
  double *corrFactor_nion_S99 = create_corrFactor_nion(endSnap, times, t0, tbreak, exponent, nion);
  
  return corrFactor_nion_S99;
}

double *get_corrFactor_nion_BPASS(int endSnap, float *times)
{
  float MyrInSec = 3.1536e13;
  float t0 = 2. * MyrInSec;
  float tbreak = 3.16 * MyrInSec;
  float exponent = -2.28;
  double nion = 9.091836e47;
  double *corrFactor_nion_BPASS = create_corrFactor_nion(endSnap, times, t0, tbreak, exponent, nion);
  
  return corrFactor_nion_BPASS;
}

double *create_corrFactor_nion(int endSnap, float *times, float t0, float tbreak, float exponent, double nion)
{
  double *corrFactor_nion = allocate_array_double(endSnap*(endSnap+1)/2, "corrFactor_nion");
  
  for(int snap=1; snap<endSnap; snap++)
  {
    for(int prevSnap=0; prevSnap<snap; prevSnap++)
    {
      corrFactor_nion[snap*(snap-1)/2 + prevSnap] = nion * calc_corrFactor_nion(times[snap]-times[prevSnap+1], times[snap]-times[prevSnap], t0, tbreak, exponent);
    }
  }
  
  return corrFactor_nion;
}

float calc_corrFactor_nion(float tinitial, float tfinal, float t0, float tbreak, float exponent)
{
  float deltat = tfinal - tinitial;
  float nion = 0.;
  
  if(tinitial >= tbreak)
  {
    nion = (pow(tfinal/t0, exponent) * tfinal / deltat - tinitial / deltat * pow(tinitial/t0, exponent)) / (1.+exponent);
  }
  else if(tfinal > tbreak)
  {
    nion = (pow(tfinal/t0, exponent) + exponent * pow(tbreak/t0, exponent) * tbreak/deltat) / (1.+exponent) - pow(tbreak/t0, exponent) * tinitial/deltat;
  }
  else
  {
    nion = pow(tbreak/t0, exponent);
  }
  
  return nion;
}