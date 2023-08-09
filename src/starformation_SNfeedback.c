#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "const.h"
#include "utils.h"

#include "dconfObj.h"
#include "gal_gtree.h"
#include "starformation_SNfeedback.h"

/*---------------------------------------------------------------------*/
/* SUPERNOVA FEEDBACK                                                  */
/*---------------------------------------------------------------------*/

float calc_SNejection_efficiency(gal_t *thisGal, dconfObj_t simParam, float corrFactorTimeStep)
{
  float G = 4.302e-6; // kpc Msol^-1 (km/s)^2
  float vc = sqrt(G * thisGal->Mvir / (thisGal->Rvir * thisGal->scalefactor));   // in km s^-1
//   printf("Mvir = %e\t rvir = %e\t vc = %e\n", thisGal->Mvir, thisGal->Rvir, vc);
  float fw = simParam->FW / (1. + (corrFactorTimeStep - 1.) * thisGal->fracMgasMer);
  float SNejection_efficiency = 0.;
  float Mgas = thisGal->MgasIni;
  
  if(simParam->delayedSNfeedback == 1)
  {
    float SNenergyPast = calc_SNenergy_past(thisGal, simParam->SNenergy);
    float factorSNenergyPast = 0.;
    if(Mgas > 0.)
      factorSNenergyPast = 1. - fw * SNenergyPast / (Mgas * vc * vc);
    
    float SNenergyCurrent = get_SNenergy_current(thisGal, simParam->SNenergy);
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

float calc_remaining_gasfraction(gal_t *thisGal, dconfObj_t simParam, float corrFactorTimeStep)
{
  float G = 4.302e-6; // kpc Msol^-1 (km/s)^2
  float vc = sqrt(G * thisGal->Mvir / (thisGal->Rvir * thisGal->scalefactor));   // in km s^-1
  float fw = simParam->FW / (1. + (corrFactorTimeStep - 1.) * thisGal->fracMgasMer);
  float remainingGasfraction = 0.;
  
  if(simParam->delayedSNfeedback == 1)
  {
    float SNenergyCurrent = get_SNenergy_current(thisGal, simParam->SNenergy);
    
    remainingGasfraction = (thisGal->fej - thisGal->feff) * (vc * vc + fw * SNenergyCurrent) / (vc * vc * (1. - thisGal->feff));
  }
  else
  {
    remainingGasfraction = thisGal->feff * (1. - thisGal->fej) / (thisGal->fej * (1. - thisGal->feff));
  }
  
  return remainingGasfraction;
}

float get_SNenergy_current(gal_t *thisGal, float *SNenergy)
{
  float result = 0.;
  int snap = thisGal->snapnumber;

  if(snap > 0)
  {
    result = SNenergy[snap*(snap-1)/2 + snap - 1];
//     result = SNenergy[snap*(snap-1)/2 + snap];
  }

  return result;
}

float calc_SNenergy_past(gal_t *thisGal, float *SNenergy)
{
  float result = 0.;
  int snap = thisGal->snapnumber;
  
  for(int prevSnap=1; prevSnap<snap; prevSnap++)
  {
    result += SNenergy[snap*(snap-1)/2 + prevSnap - 1] * thisGal->stellarmasshistory[prevSnap];
//     result += SNenergy[snap*(snap-1)/2 + prevSnap] * thisGal->stellarmasshistory[prevSnap];
  }
  
  return result;
}

/*---------------------------------------------------------------------*/
/* COMPUTING DELAYED SUPERNOVA FEEDBACK TABLE                          */
/*---------------------------------------------------------------------*/

float *get_SNenergy_delayed(dconfObj_t simParam)
{
  int endSnap = simParam->endSnap;
  float *times = simParam->timeSnaps;
  int burstySF = simParam->burstySF;
  float *SNenergy_delayed = allocate_array_float(endSnap*(endSnap+1)/2, "SNenergy_delayed");
  double invHubble_h = 1./simParam->hubble_h;
  
  for(int snap=1; snap<endSnap; snap++)
  {
    for(int prevSnap=0; prevSnap<snap; prevSnap++)
    {
      SNenergy_delayed[snap*(snap-1)/2 + prevSnap] = 5.e7 * invHubble_h * get_SNfraction_per_Msun(times[snap], times[snap-1], times[prevSnap+1], times[prevSnap], burstySF);   // factor due to conversion form erg to Msun (km/s)^2
    }
  }
  
  return SNenergy_delayed;
}

float get_SNfraction_per_Msun(float tsnap, float tsnapPrev, float tprevSnap, float tprevprevSnap, int burstySF)
{
  float IMFslope = -2.35;
  float IMFminMstar = 0.1;
  float IMFmaxMstar = 100.;
  float SNminMstar = 8.;
  
  float deltaTprevSnap = tsnap - tprevSnap;
  float deltaTprevprevSnap = tsnap - tprevprevSnap;
    
  float MprevSnap = IMFmaxMstar;
  if(tsnap != tprevSnap)
    MprevSnap = get_SNmass(deltaTprevSnap*secInMyr);
  float MprevprevSnap = get_SNmass(deltaTprevprevSnap*secInMyr);

  if(MprevSnap > IMFmaxMstar)
    MprevSnap = IMFmaxMstar;
  if(MprevprevSnap > IMFmaxMstar)
    MprevprevSnap = IMFmaxMstar;
  if(MprevSnap < SNminMstar)
    MprevSnap = SNminMstar;
  if(MprevprevSnap < SNminMstar)
    MprevprevSnap = SNminMstar;
  
  float SNfraction_per_Msun = 0.;
//   float SNfraction_per_Msun_burstySF = 0.;
  
  if(burstySF == 1)
  {
    if(MprevprevSnap < MprevSnap)
      SNfraction_per_Msun = calc_SNfraction_per_Msun(MprevprevSnap, MprevSnap, IMFminMstar, IMFmaxMstar, IMFslope);
  }
  else
  {
//     if(MprevprevSnap < MprevSnap)
//       SNfraction_per_Msun_burstySF = calc_SNfraction_per_Msun(MprevprevSnap, MprevSnap, IMFminMstar, IMFmaxMstar, IMFslope);
    SNfraction_per_Msun = calc_SNfraction_per_Msun_contSF(tsnap, tsnapPrev, tprevSnap, tprevprevSnap, SNminMstar, IMFminMstar, IMFmaxMstar, IMFslope);
  }
  
  return SNfraction_per_Msun;
}

float get_SNmass(float timeInMyr)
{
  float SNmass = pow(0.83333e-3*(timeInMyr - 3.), -0.540541);
  
  return SNmass;
}

float get_SNtime(float massInMsun)
{
  float SNtime = 1.2e3 * pow(massInMsun, -1.85) + 3.;
  
  return SNtime;
}

float calc_SNfraction_per_Msun(float MSNlow, float MSNup, float Mlow, float Mup, float slope)
{
  float SNfraction_per_Msun = (2.+slope)/(1.+slope) * (pow(MSNup, 1.+slope) - pow(MSNlow, 1.+slope)) / (pow(Mup, 2.+slope) - pow(Mlow, 2.+slope));
  
  return SNfraction_per_Msun;
}

float calc_SNfraction_per_Msun_contSF(float tsnap, float tsnapPrev, float tprevSnap, float tprevprevSnap, float MSNlow, float Mlow, float Mup, float slope)
{
  float tSNlow = get_SNtime(MSNlow);
  float tup = get_SNtime(Mup);
  
  float exponent = 0.540541 * (-slope-1.);
  float preFactor = (2.+slope) / (1.+slope) * pow(0.83333e-3, exponent) / (pow(Mlow, 2.+slope) - pow(Mup, 2.+slope));
  float SFR = 1. / (tprevSnap - tprevprevSnap) * MyrInSec; // in Msun/sec
  
  /* tmin part of the integration */   
  float resultTminLow = 0.;
  float tminLowLow = MAX(tsnapPrev*secInMyr, tprevprevSnap*secInMyr + tup);
  float tminLowUp = MIN(tsnap*secInMyr, tprevprevSnap*secInMyr + tSNlow); 
  if(tminLowLow < tminLowUp)
    resultTminLow = (pow(tminLowUp - tprevprevSnap*secInMyr - 3., exponent + 1.) - pow(tminLowLow - tprevprevSnap*secInMyr - 3., exponent + 1.)) / (exponent + 1.);
  
  float resultTminUp = 0.;
  float tminUpLow = MAX(tsnapPrev*secInMyr, tprevprevSnap*secInMyr + tSNlow);
  float tminUpUp = MIN(tsnap*secInMyr, tprevSnap*secInMyr + tSNlow);
  if(tminUpLow < tminUpUp)
    resultTminUp = pow(tSNlow - 3., exponent) * (tminUpUp - tminUpLow);
  
  /* tmax part of the integration */   
  float resultTmaxLow = 0.;
  float tmaxLowLow = MAX(tsnapPrev*secInMyr, tprevprevSnap*secInMyr + tup);
  float tmaxLowUp = MIN(tsnap*secInMyr, tprevSnap*secInMyr + tup);
  if(tmaxLowLow < tmaxLowUp)
    resultTmaxLow = pow(tup - 3., exponent) * (tmaxLowUp - tmaxLowLow);
  
  float resultTmaxUp = 0.;
  float tmaxUpLow = MAX(tsnapPrev*secInMyr, tprevSnap*secInMyr + tup);
  float tmaxUpUp = MIN(tsnap*secInMyr, tprevSnap*secInMyr + tSNlow);
  if(tmaxUpLow < tmaxUpUp)
    resultTmaxUp = (pow(tmaxUpUp - tprevSnap*secInMyr - 3., exponent + 1.) - pow(tmaxUpLow - tprevSnap*secInMyr - 3., exponent + 1.)) / (exponent + 1.);
  
  float result = SFR * preFactor * (resultTminLow + resultTminUp - resultTmaxLow - resultTmaxUp);
    
  return result;
}
