#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include "dconfObj.h"
#include "gal_gtree.h"
#include "starformation_SNfeedback.h"
#include "fesc.h"

double get_fesc(gal_t *thisGal, dconfObj_t simParam)
{
  double fesc = 1.;
  int fesc_model = simParam->fesc_model;
  
  if(fesc_model == 2)
    fesc = get_fesc_MHinc(thisGal, simParam);
  else if(fesc_model == 1)
    fesc = get_fesc_MHdec(thisGal, simParam);
  else if(fesc_model == 3)
    fesc = get_fesc_SN(thisGal, simParam);
  else 
    fesc = get_fesc_const(simParam);
  
  return fesc;
}

double get_fesc_const(dconfObj_t simParam)
{
  double fesc = simParam->fesc;
  
  return fesc;
}

double get_fesc_MHinc(gal_t *thisGal, dconfObj_t simParam)
{
  double Mvir = thisGal->Mvir;
  double MvirLow = simParam->MHlow;
  double MvirHigh = simParam->MHhigh;
  double fescLow = simParam->fescLow;
  double fescHigh = simParam->fescHigh;
  
  double invMvirLow = 1./MvirLow;
  
  double fesc = fescLow * pow(fescLow/fescHigh, - log10(Mvir*invMvirLow) / log10(MvirHigh*invMvirLow));
  
  return fesc;
}

double get_fesc_MHdec(gal_t *thisGal, dconfObj_t simParam)
{
  double Mvir = thisGal->Mvir;
  double MvirLow = simParam->MHlow;
  double MvirHigh = simParam->MHhigh;
  double fescLow = simParam->fescLow;
  double fescHigh = simParam->fescHigh;
  
  double invMvirLow = 1./MvirLow;
  
  double fesc = (1.-fescLow) * pow((1.-fescLow)/(1.-fescHigh), - log10(Mvir*invMvirLow) / log10(MvirHigh*invMvirLow));
  
  return fesc;
}

double get_fesc_SN(gal_t *thisGal, dconfObj_t simParam)
{
  double fej = calc_SNejection_efficiency(thisGal, simParam, thisGal->MgasIni);
  double feff = thisGal->feff;
  
  double fesc = 1.;
  if(fej > 0.)
    fesc = feff/fej;
  if(fesc > 1.)
    fesc = 1.;
  fesc = fesc * simParam->fesc;
  
  return fesc;
}