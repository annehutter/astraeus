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
#include "timestep_model.h"
#include "evol_gal.h"

/*---------------------------------------------------------------------*/
/* This function evolves the properties of a galaxy for one time step  */
/*---------------------------------------------------------------------*/
void evolve_gal(gal_t *thisGal, gtree_t *thisGtree, dconfObj_t simParam)
{
  int32_t descGalID = thisGal->localDescID;
  if(descGalID == -1)
        descGalID = 0;
  gal_t *descGal = &(thisGtree->galaxies[descGalID]);
      
  /* adjust DM mass to continuously rise */
  correct_Mvir_to_continously_rise(thisGal);
  
  /* get gas from inflows & progenitors */
  get_accreted_and_merged_gas(simParam, thisGal);
  
  /* radiative feedback (how much gas can the galaxy keep) */
  apply_radFeedback_and_update_zreion(simParam, thisGal, descGal);
  
  /* saving initial gas mass (before SN feedback) */
  store_available_gas(thisGal);
  
  /* star formation & SN feedback */
  do_starFormation_and_SNfeedback(simParam, thisGal, descGal);

  /* copying properties to be inherited to descendant */
  copy_properties_to_descGal(simParam, thisGal, descGal);
}

/*---------------------------------------------------------------------*/
/* Initialising the stellar mass history of a galaxy                   */
/*---------------------------------------------------------------------*/
void start_stellarmasshistory(gal_t *thisGal)
{
  if(thisGal->stellarmasshistory == NULL)
    thisGal->stellarmasshistory = allocate_array_float(thisGal->snapnumber+1, "stellarmasshistory");
}

/*---------------------------------------------------------------------*/
/* Storing the stellar mass history of a galaxy and copying it to the  */
/* galaxy's descendant                                                 */
/*---------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------*/
/* This function frees the stellar mass history of a galaxy            */
/*---------------------------------------------------------------------*/
void clean_gal(gal_t *thisGal, int outSnap)
{
  if(thisGal->snapnumber != outSnap)
    free(thisGal->stellarmasshistory);
}

/*---------------------------------------------------------------------*/
/* Analytic derivation for the virial radius from the virial mass      */
/*---------------------------------------------------------------------*/
float calc_Rvir(float Mvir, double omega_m)
{
  float result = 0.;
  
  result = 0.0169 * pow(Mvir/omega_m, 0.333333);
  
  return result;
}

/*---------------------------------------------------------------------*/
/* Correct Mvir & Rvir if Mvir is not continuously rising in trees     */
/*---------------------------------------------------------------------*/
void correct_Mvir_to_continously_rise(gal_t *thisGal)
{
  //   thisGal->Rvir = calc_Rvir(thisGal->Mvir, simParam->omega_m);

  if(thisGal->Mvir < thisGal->Mvir_prog)
  {
    thisGal->Rvir = thisGal->Rvir * pow(thisGal->Mvir_prog/thisGal->Mvir, 0.3333);
    thisGal->Mvir = thisGal->Mvir_prog;
  }
}

/*---------------------------------------------------------------------*/
/* Get accreted & merged gas mass & fraction of merged gas             */
/*---------------------------------------------------------------------*/
void get_accreted_and_merged_gas(dconfObj_t simParam, gal_t *thisGal)
{
  float MgasMer = 0.;
  if(thisGal->numProg == 0)
  {
    thisGal->Mgas = thisGal->Mvir * simParam->omega_b / simParam->omega_m;
    thisGal->fracMgasMer = 0.;
  }
  else if(thisGal->Mvir_prog < thisGal->Mvir)
  {
    MgasMer = thisGal->Mgas;
    thisGal->Mgas += (thisGal->Mvir - thisGal->Mvir_prog) * simParam->omega_b / simParam->omega_m;
    thisGal->fracMgasMer = MgasMer / thisGal->Mgas;
  }
}

/*---------------------------------------------------------------------*/
/* Compute reduction of gas mass for star formation due to radiative   */
/* feedback and track when a galaxy's environment was reionized        */
/*---------------------------------------------------------------------*/
void apply_radFeedback_and_update_zreion(dconfObj_t simParam, gal_t *thisGal, gal_t * descGal)
{
  float fg = 1.;
  if(simParam->radfeedback == 1)
  {
    fg = fmin(1., calc_radiative_feedback(thisGal, descGal, simParam));
    thisGal->fg = fg;
    if(fg > 1.) printf("fg = %e\n", fg);
    if(thisGal->Mgas > fg * simParam->omega_b / simParam->omega_m * thisGal->Mvir)
        thisGal->Mgas = fg * simParam->omega_b / simParam->omega_m * thisGal->Mvir;
  }
  else if(simParam->reion == 1)
  {
    pass_zreion_to_descGal(thisGal, descGal, simParam);
  }
}

/*---------------------------------------------------------------------*/
/* Store gas mass available for star formation                         */
/*---------------------------------------------------------------------*/
void store_available_gas(gal_t *thisGal)
{
  thisGal->MgasIni = thisGal->Mgas;
}

/*---------------------------------------------------------------------*/
/* Compute stellar mass formed, track stellar mass history and gas     */
/* ejected due to SN feedback                                          */
/*---------------------------------------------------------------------*/
void do_starFormation_and_SNfeedback(dconfObj_t simParam, gal_t *thisGal, gal_t *descGal)
{
  float fej = 1.;
  float corrFactorTimeStep = get_corrFactorTimeStep(simParam, thisGal->snapnumber);
  float fstar = simParam->FS * corrFactorTimeStep;
  float feff = fstar;
  float newMstar = 0.;
  
  /* supernova feedback */
  start_stellarmasshistory(thisGal);
  fej = calc_SNejection_efficiency(thisGal, simParam, corrFactorTimeStep);
 
  /* star formation */
  feff = fmin(fstar, fej);
  newMstar = thisGal->Mgas*feff;
  thisGal->feff = feff;

  /* total stellar mass and final gas mass */
  thisGal->Mstar += newMstar;
  if(fej > 0.)
    thisGal->Mgas = (thisGal->Mgas - newMstar)*(1 - feff/fej);
  else
    thisGal->Mgas = 0.;

  /* storing stellar mass history of the galaxy */
  track_stellarmasshistory(thisGal, descGal, simParam, newMstar);
}

/*---------------------------------------------------------------------*/
/* Copy stellar, gas and DM to descendant galaxy (required to compute  */
/* accreted DM and gas masses and merged gas mass)                     */
/*---------------------------------------------------------------------*/
void copy_properties_to_descGal(dconfObj_t simParam, gal_t *thisGal, gal_t *descGal)
{
  if((descGal->snapnumber != simParam->endSnap) && (thisGal->localDescID != -1))
  {
    descGal->Mvir_prog += thisGal->Mvir;
    descGal->Mgas += thisGal->Mgas;
    descGal->Mstar += thisGal->Mstar;
  }
}

/*---------------------------------------------------------------------*/
/* This function passes the redshift of reionization the galaxy's      */
/* descendant. This is required when radiative feedback is disabled.   */
/*---------------------------------------------------------------------*/
void pass_zreion_to_descGal(gal_t *thisGal, gal_t *descGal, dconfObj_t simParam)
{
  float z = 1./thisGal->scalefactor - 1.;
  float zreion = thisGal->zreion;
  
  if(simParam->reion == 1)
  {
    if((descGal->snapnumber != simParam->endSnap) && (thisGal->localDescID != -1))
    {
      if(descGal->zreion < zreion && zreion >= z)
      {
        descGal->zreion = zreion;
      }
    }
  }
}
