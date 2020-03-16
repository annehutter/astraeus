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

#include "utils.h"
#include "dconfObj.h"
#include "gal_gtree.h"
#include "radiative_feedback.h"
#include "starformation_SNfeedback.h"
#include "evol_gal.h"

void evolve_gal(gal_t *thisGal, gtree_t *thisGtree, dconfObj_t simParam)
{
  int32_t descGalID = thisGal->localDescID;
  if(descGalID == -1)
  	descGalID = 0;
  gal_t *descGal = &(thisGtree->galaxies[descGalID]);
  float fg = 1.;
  float fej = 1.;
  float fstar = simParam->FS;
  float feff = fstar;
  float newMstar = 0.;
     
  /* adjust DM mass to continuously rise */
  if(thisGal->Mvir < thisGal->Mvir_prog)
  {
    thisGal->Rvir = thisGal->Rvir * pow(thisGal->Mvir_prog/thisGal->Mvir, 0.3333);
    thisGal->Mvir = thisGal->Mvir_prog;
  }
  
  /* get gas from inflows & progenitors */
  if(thisGal->numProg == 0)
    thisGal->Mgas = thisGal->Mvir * simParam->omega_b / simParam->omega_m;
  else if(thisGal->Mvir_prog < thisGal->Mvir)
    thisGal->Mgas += (thisGal->Mvir - thisGal->Mvir_prog) * simParam->omega_b / simParam->omega_m;
  
  /* radiative feedback (how much gas can the galaxy keep) */
  if(simParam->radfeedback == 1)
  {
    fg = fmin(1., calc_radiative_feedback(thisGal, descGal, simParam));
    thisGal->fg = fg;
    if(fg > 1.) printf("fg = %e\n", fg);
    if(thisGal->Mgas > fg * simParam->omega_b / simParam->omega_m * thisGal->Mvir)
        thisGal->Mgas = fg * simParam->omega_b / simParam->omega_m * thisGal->Mvir;
  }
  
  thisGal->MgasIni = thisGal->Mgas;
  
  /* Supernova */
  start_stellarmasshistory(thisGal);
  fej = calc_SNejection_efficiency(thisGal, simParam, thisGal->Mgas);
  
  /* star formation */
  feff = fmin(fstar, fej);
  newMstar = thisGal->Mgas*feff;
  thisGal->feff = feff;

  /* book keeping */
  thisGal->Mstar += newMstar;
  if(fej > 0.)
    thisGal->Mgas = (thisGal->Mgas - newMstar)*(1 - feff/fej);
  else
    thisGal->Mgas = 0.;
  
  track_stellarmasshistory(thisGal, descGal, simParam, newMstar);

  if((descGal->snapnumber != simParam->endSnap) && (thisGal->localDescID != -1))
  {
    descGal->Mvir_prog += thisGal->Mvir;
    descGal->Mgas += thisGal->Mgas;
    descGal->Mstar += thisGal->Mstar;
  }
}

void start_stellarmasshistory(gal_t *thisGal)
{
  if(thisGal->stellarmasshistory == NULL)
    thisGal->stellarmasshistory = allocate_array_float(thisGal->snapnumber+1, "stellarmasshistory");
}

void track_stellarmasshistory(gal_t *thisGal, gal_t *descGal, dconfObj_t simParam, float newMstar)
{
  if(thisGal->stellarmasshistory == NULL)
    thisGal->stellarmasshistory = allocate_array_float(thisGal->snapnumber+1, "stellarmasshistory");
  thisGal->stellarmasshistory[thisGal->snapnumber] = newMstar;

  if((descGal->snapnumber != simParam->endSnap) && (thisGal->localDescID != -1))
  {
    if(descGal->stellarmasshistory == NULL)
      descGal->stellarmasshistory = allocate_array_float(descGal->snapnumber+1, "stellarmasshistory");

    for(int i=0; i<thisGal->snapnumber+1; i++)
      descGal->stellarmasshistory[i] += thisGal->stellarmasshistory[i];
  }
}

void clean_gal(gal_t *thisGal, int outSnap)
{
  if(thisGal->snapnumber != outSnap)
    free(thisGal->stellarmasshistory);
}
