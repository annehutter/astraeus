#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"

#include "dconfObj.h"
#include "SNfeedback.h"

/*---------------------------------------------------------------------*/
/* SUPERNOVA FEEDBACK                                                  */
/*---------------------------------------------------------------------*/

double *get_SNejection_efficiency(dconfObj_t simParam, int numGal, double *MvirDivRvir, double *MgasIni, double *stellarMassHistory, int currSnap)
{
  double *SNejectionEfficiency = allocate_array_double(numGal, "SNejectionEfficiency");
  
  for(int gal=0; gal<numGal; gal++)
  {
    SNejectionEfficiency[gal] = calc_SNejection_efficiency(simParam, MvirDivRvir[gal], MgasIni[gal], &(stellarMassHistory[gal*(currSnap + 1)]), currSnap);
  }
  
  return SNejectionEfficiency;
}

double calc_SNejection_efficiency(dconfObj_t simParam, double MvirDivRvir, double Mgas, double *stellarmasshistory, int currSnap)
{
  float scalefactor = simParam->scalefactors[currSnap];
  float G = 4.302e-6; // kpc Msol^-1 (km/s)^2
  float vc = sqrt(G * MvirDivRvir / scalefactor);   // in km s^-1
  float fw = simParam->fw;
  double SNejection_efficiency = 0.;
  
  if(simParam->delayedSNfeedback == 1)
  {
    float SNenergyPast = calc_SNenergy_past(currSnap, simParam->SNenergy, stellarmasshistory);
    float factorSNenergyPast = 0.;
    if(Mgas > 0.)
      factorSNenergyPast = 1. - fw * SNenergyPast / (Mgas * vc * vc);
    
    float SNenergyCurrent = get_SNenergy_current(currSnap, simParam->SNenergy);
    float factorSNenergyCurrent = vc * vc / (vc * vc + fw * SNenergyCurrent);
    
    SNejection_efficiency = factorSNenergyCurrent * fmax(0., factorSNenergyPast);
  }
  else
  {
    float vs = 611.; // km/s
    SNejection_efficiency = vc * vc / (vc * vc + fw * vs * vs);
  }
  
  return SNejection_efficiency;
}

float get_SNenergy_current(int currSnap, float *SNenergy)
{
  float result = 0.;
  int snap = currSnap;

  if(snap > 0)
    result = SNenergy[snap*(snap-1)/2 + snap - 1];
  
  return result;
}

float calc_SNenergy_past(int currSnap, float *SNenergy, double *stellarmasshistory)
{
  float result = 0.;
  int snap = currSnap;
  
  for(int prevSnap=1; prevSnap<snap; prevSnap++)
  {
    result += SNenergy[snap*(snap-1)/2 + prevSnap - 1] * stellarmasshistory[prevSnap];
  }
  
  return result;
}

/*---------------------------------------------------------------------*/
/* COMPUTING DELAYED SUPERNOVA FEEDBACK TABLE                          */
/*---------------------------------------------------------------------*/

float *get_SNenergy_delayed(int endSnap, float *times)
{
  float *SNenergy_delayed = allocate_array_float(endSnap*(endSnap+1)/2, "SNenergy_delayed");
  
  for(int snap=0; snap<endSnap; snap++)
  {
    for(int prevSnap=0; prevSnap<snap; prevSnap++)
    {
      SNenergy_delayed[snap*(snap-1)/2 + prevSnap] = 5.e7 * get_SNfraction_per_Msun(times[snap], times[prevSnap+1], times[prevSnap]);   // factor due to conversion form erg to Msun (km/s)^2
    }
  }
  
  return SNenergy_delayed;
}

float get_SNfraction_per_Msun(float tsnap, float tprevSnap, float tprevprevSnap)
{
  float IMFslope = -2.35;
  float IMFminMstar = 0.1;
  float IMFmaxMstar = 100.;
  float SNminMstar = 8.;
  
  float deltaTprevSnap = tsnap - tprevSnap;
  float deltaTprevprevSnap = tsnap - tprevprevSnap;
    
  float secInGyr = 3.170979e-17;
  float MprevSnap = IMFmaxMstar;
  if(tsnap != tprevSnap)
    MprevSnap = get_SNmass(deltaTprevSnap*secInGyr);
  float MprevprevSnap = get_SNmass(deltaTprevprevSnap*secInGyr);

  if(MprevSnap > IMFmaxMstar)
    MprevSnap = IMFmaxMstar;
  if(MprevprevSnap > IMFmaxMstar)
    MprevprevSnap = IMFmaxMstar;
  if(MprevSnap < SNminMstar)
    MprevSnap = SNminMstar;
  if(MprevprevSnap < SNminMstar)
    MprevprevSnap = SNminMstar;
  
  float SNfraction_per_Msun = 0.;
  if(MprevprevSnap < MprevSnap)
    SNfraction_per_Msun = calc_SNfraction_per_Msun(MprevprevSnap, MprevSnap, IMFminMstar, IMFmaxMstar, IMFslope);
  
  return SNfraction_per_Msun;
}

float get_SNmass(float timeInGyr)
{
  float SNmass = pow(0.83333*(timeInGyr - 0.003), -0.540541);
  
  return SNmass;
}

float calc_SNfraction_per_Msun(float MSNlow, float MSNup, float Mlow, float Mup, float slope)
{
  float SNfraction_per_Msun = (2.+slope)/(1.+slope) * (pow(MSNup, 1.+slope) - pow(MSNlow, 1.+slope)) / (pow(Mup, 2.+slope) - pow(Mlow, 2.+slope));
  
  return SNfraction_per_Msun;
}