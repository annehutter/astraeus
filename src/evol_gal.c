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
#include "radiative_feedback.h"
#include "starformation_SNfeedback.h"
#include "timestep_model.h"
#include "evol_gal.h"
#include "metallicity.h"

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
//   printf("Mvir = %e\t Rvir = %e\t Rvir(recalc) = %e\n", thisGal->Mvir, thisGal->Rvir, calc_Rvir(thisGal->Mvir, simParam->omega_m));
  correct_Mvir_to_continously_rise(thisGal);
  
  /* get metals from inflows & progenitors */
#if defined WITHMETALS
  if(simParam->metals == 1) 
    get_accreted_and_merged_metals(simParam, thisGal);
#endif

  /* get gas from inflows & progenitors */
  get_accreted_and_merged_gas(simParam, thisGal);
  
  /* radiative feedback (how much gas can the galaxy keep) */
  apply_radFeedback_and_update_zreion(simParam, thisGal, descGal);

  /* apply gas loss due to radiative feedback to metals */
#if defined WITHMETALS
  if(simParam->metals == 1)
    apply_radFeedback_to_metals(simParam, thisGal);
#endif
  
  /* saving initial gas mass (before SN feedback) */
  store_available_gas(thisGal);
  
  /* saving initial metal masses (before SN feedback) */
#if defined WITHMETALS
  if(simParam->metals == 1)
    store_available_metals(thisGal);
#endif

  /* star formation & SN feedback */
  do_starFormation_and_SNfeedback(simParam, thisGal, descGal);

#if defined WITHMETALS
  if(simParam->metals == 1)
    produce_and_eject_metals(simParam, thisGal, descGal);
#endif

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

#if defined WITHMETALS
void start_metalmasshistory(gal_t *thisGal)
{
  if(thisGal->metalmasshistory == NULL)
    thisGal->metalmasshistory = allocate_array_float(thisGal->snapnumber+1, "metalmasshistory");
}
#endif

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

#if defined WITHMETALS
void track_metalmasshistory(gal_t *thisGal, gal_t *descGal, dconfObj_t simParam, float newMmetal)
{
  if(thisGal->metalmasshistory == NULL)
    thisGal->metalmasshistory = allocate_array_float(thisGal->snapnumber+1, "metalmasshistory");
  thisGal->metalmasshistory[thisGal->snapnumber] = newMmetal;

  if((descGal->snapnumber != simParam->endSnap) && (thisGal->localDescID != -1))
  {
    if(descGal->metalmasshistory == NULL)
      descGal->metalmasshistory = allocate_array_float(descGal->snapnumber+1, "metalmasshistory");

    for(int i=0; i<thisGal->snapnumber+1; i++)
      descGal->metalmasshistory[i] += thisGal->metalmasshistory[i];
  }
}
#endif

/*---------------------------------------------------------------------*/
/* This function frees the stellar mass history of a galaxy            */
/*---------------------------------------------------------------------*/
void clean_gal(dconfObj_t simParam, gal_t *thisGal, int outSnap)
{
  if(thisGal->snapnumber != outSnap)
  {
    free(thisGal->stellarmasshistory);
#if defined WITHMETALS
    if(simParam->metals == 1)
      free(thisGal->metalmasshistory);
#endif
  }
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
/* Get accreted & merged metal mass                                    */
/*---------------------------------------------------------------------*/
#if defined WITHMETALS
void get_accreted_and_merged_metals(dconfObj_t simParam, gal_t *thisGal)
{
  float MgasAcc = 0.;
  float MmetalAcc = 0., MmetalAccO = 0., MmetalAccFe = 0.;
  float metallicityIgm = 0., metallicityIgmO = 0., metallicityIgmFe = 0.;
  float dustfractionIgm = 0.;

  /* Compute average IGM metallicity from previous snap */
  if(thisGal->snapnumber != 0)
  {
    metallicityIgm   = thisGal->igmMetallicity[0]; //simParam->meanMetallicityIgm[thisGal->snapnumber-1];
    metallicityIgmO  = thisGal->igmMetallicity[1]; //simParam->meanMetallicityIgmO[thisGal->snapnumber-1];
    metallicityIgmFe = thisGal->igmMetallicity[2]; //simParam->meanMetallicityIgmFe[thisGal->snapnumber-1];
    dustfractionIgm = thisGal->igmDustFraction; //simParam->meanDustfractionIgm[thisGal->snapnumber-1];
  }

  /* Get gas from INFLOWS and PROGENITORS */
  if(thisGal->numProg == 0)
  {
    MgasAcc            = thisGal->Mvir * simParam->omega_b / simParam->omega_m;
    MmetalAcc         = metallicityIgm   * MgasAcc;
    MmetalAccO        = metallicityIgmO  * MgasAcc;
    MmetalAccFe       = metallicityIgmFe * MgasAcc;
    thisGal->Mmetal[0] = MmetalAcc;
    thisGal->Mmetal[1] = MmetalAccO;
    thisGal->Mmetal[2] = MmetalAccFe;
    thisGal->Mdust  = dustfractionIgm   * MgasAcc;
  }
  else if(thisGal->Mvir_prog < thisGal->Mvir)
  {
    MgasAcc           = (thisGal->Mvir - thisGal->Mvir_prog) * simParam->omega_b / simParam->omega_m;
    MmetalAcc           = metallicityIgm   * MgasAcc;
    MmetalAccO          = metallicityIgmO  * MgasAcc;
    MmetalAccFe         = metallicityIgmFe * MgasAcc;
    thisGal->Mmetal[0] += MmetalAcc;
    thisGal->Mmetal[1] += MmetalAccO;
    thisGal->Mmetal[2] += MmetalAccFe;
    thisGal->Mdust  += dustfractionIgm   * MgasAcc;
  }
  
  thisGal->fracMmetalMer[0] = 1. - (MmetalAcc    / thisGal->Mmetal[0]);
  thisGal->fracMmetalMer[1] = 1. - (MmetalAccO  / thisGal->Mmetal[1]);
  thisGal->fracMmetalMer[2] = 1. - (MmetalAccFe / thisGal->Mmetal[2]);
}
#endif

/*---------------------------------------------------------------------*/
/* Compute reduction of gas mass for star formation due to radiative   */
/* feedback and track when a galaxy's environment was reionized        */
/*---------------------------------------------------------------------*/
void apply_radFeedback_and_update_zreion(dconfObj_t simParam, gal_t *thisGal, gal_t * descGal)
{
  float fg = 1.;
  float MgasMax = 0.;
  
  if(simParam->radfeedback == 1)
  {
    fg = fmin(1., calc_radiative_feedback(thisGal, descGal, simParam));
    thisGal->fg = fg;
    assert(fg <= 1.);
    MgasMax = thisGal->fg * simParam->omega_b / simParam->omega_m * thisGal->Mvir;
    if(thisGal->Mgas > MgasMax)
      thisGal->Mgas = MgasMax;
  }
  else if(simParam->reion == 1)
  {
    pass_zreion_to_descGal(thisGal, descGal, simParam);
  }
}

#if defined WITHMETALS
void apply_radFeedback_to_metals(dconfObj_t simParam, gal_t *thisGal)
{
  float metallicityIni = 0., metallicityIniO = 0., metallicityIniFe = 0.;
  float dustfractionIni = 0.;
  float MgasMax = 0.;
  
  if(simParam->radfeedback == 1)
  {
    if(thisGal->Mgas > 0.) 
    {
      metallicityIni    = thisGal->Mmetal[0] / thisGal->Mgas;
      metallicityIniO  = thisGal->Mmetal[1] / thisGal->Mgas;
      metallicityIniFe = thisGal->Mmetal[2] / thisGal->Mgas;
      dustfractionIni    = thisGal->Mdust  / thisGal->Mgas; 
    }
    
    MgasMax = thisGal->fg * simParam->omega_b / simParam->omega_m * thisGal->Mvir;
    if(thisGal->Mgas > MgasMax)
    {
        thisGal->Mmetal[0] = metallicityIni * MgasMax;
        thisGal->Mmetal[1] = metallicityIniO * MgasMax;
        thisGal->Mmetal[2] = metallicityIniFe * MgasMax; 
        thisGal->Mdust  = dustfractionIni * MgasMax;
    }
  }
}
#endif

/*---------------------------------------------------------------------*/
/* Store gas mass available for star formation                         */
/*---------------------------------------------------------------------*/
void store_available_gas(gal_t *thisGal)
{
  thisGal->MgasIni = thisGal->Mgas;
}

#if defined WITHMETALS
void store_available_metals(gal_t *thisGal)
{
  thisGal->MmetalIni[0] = thisGal->Mmetal[0];
  thisGal->MmetalIni[1] = thisGal->Mmetal[1];
  thisGal->MmetalIni[2] = thisGal->Mmetal[2];
}
#endif

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
  thisGal->fej = fej;
 
  /* star formation */
  feff = fmin(fstar, fej);
  newMstar = thisGal->Mgas*feff;
  thisGal->feff = feff;

  /* total stellar mass and final gas mass */
  thisGal->Mstar += newMstar;
  if(fej > 0.)
//     thisGal->Mgas = (thisGal->Mgas - newMstar)*(1 - feff/fej);
    thisGal->Mgas = (thisGal->Mgas - newMstar) * calc_remaining_gasfraction(thisGal, simParam, corrFactorTimeStep);
  else
    thisGal->Mgas = 0.;

  /* storing stellar mass history of the galaxy */
  track_stellarmasshistory(thisGal, descGal, simParam, newMstar);
}

#if defined WITHMETALS
void produce_and_eject_metals(dconfObj_t simParam, gal_t *thisGal, gal_t *descGal)
{
  int snap = thisGal->snapnumber;
  float metallicityIni = 0., metallicityIniO = 0., metallicityIniFe = 0., dustfractionIni = 0.;
  float metallicity = 0., metallicityO = 0., metallicityFe = 0., dustfraction = 0.;
  float MmetalEj = 0., MmetalEjO = 0., MmetalEjFe = 0., MdustEj = 0.;
  float newMstar = thisGal->stellarmasshistory[snap];
  float sfr = 0.;
  float tmpMgas = thisGal->MgasIni, MgasEj = 0.;

  float G_t       = 0.;    /* (G(t) x Deltat) in Ucci et al. (2021) */
  float e_Z       = 0.;    /* (e_Z  x Deltat) in Ucci et al. (2021) */
  float e_Z_O     = 0.;    /* (e_Z  x Deltat) in Ucci et al. (2021) for OXYGEN */
  float e_Z_Fe    = 0.;    /* (e_Z  x Deltat) in Ucci et al. (2021) for IRON */
  float y_d       = 0.;    /* (y_d  x Deltat) in Ucci et al. (2021) for dust */
  float zetaej    = 0.;    /* metal loading factor for ejection */
  
  /* Compute initial metallicities */
  if(thisGal->MgasIni != 0)
  {
    metallicityIni   = thisGal->Mmetal[0] / thisGal->MgasIni;
    metallicityIniO  = thisGal->Mmetal[1] / thisGal->MgasIni;
    metallicityIniFe = thisGal->Mmetal[2] / thisGal->MgasIni;
    dustfractionIni  = thisGal->Mdust / thisGal->MgasIni;
  }
  
  start_metalmasshistory(thisGal);
  
  /* Compute SFR in galaxy [ M_sun yr^-1 ]  */
  if(snap > 0)
    sfr = YrInSec * newMstar * simParam->invdeltat[snap];

  /* Compute metal and gas masses from SN explosions */
  tmpMgas -= newMstar; 
  thisGal->Mmetal[0] -= metallicityIni   * newMstar;
  thisGal->Mmetal[1] -= metallicityIniO  * newMstar;
  thisGal->Mmetal[2] -= metallicityIniFe * newMstar;
  thisGal->Mdust     -= dustfractionIni  * newMstar; 
  compute_metal_mass(thisGal, simParam, sfr, (1. / simParam->invdeltat[snap]), &G_t, &e_Z, &e_Z_O, &e_Z_Fe, &y_d);
  
  if(snap < 1) /* fix when gsl integration or so fails, returning NAN */
  {
    if(isnan(G_t)) G_t = 0.;
    if(isnan(e_Z)) e_Z = 0.;
    if(isnan(e_Z_O)) e_Z_O = 0.;
    if(isnan(e_Z_Fe)) e_Z_Fe = 0.;
    if(isnan(y_d)) y_d = 0.;
  }

  tmpMgas += G_t + e_Z;
  thisGal->Mmetal[0] += e_Z;
  thisGal->Mmetal[1] += e_Z_O;
  thisGal->Mmetal[2] += e_Z_Fe;
  thisGal->Mdust     += y_d;

  // Assign metal and gass masses released to outputs
  thisGal->MgasNew   = G_t + e_Z;
  thisGal->MmetalNew[0] = e_Z;
  thisGal->MmetalNew[1] = e_Z_O;
  thisGal->MmetalNew[2] = e_Z_Fe; 
  
  /* get metallicities before ejection */
  metallicity   = thisGal->Mmetal[0] / tmpMgas;
  metallicityO  = thisGal->Mmetal[1] / tmpMgas;
  metallicityFe = thisGal->Mmetal[2] / tmpMgas;
  dustfraction  = thisGal->Mdust / tmpMgas;
  
//   if(thisGal->Mvir > 5.e9) printf("snap = %d: Mdust = %e\t y_d = %e\t Mgas = %e\n", snap, thisGal->Mdust, y_d, thisGal->Mgas);
  
  /* the ejected masses */
//   float feffDIVfej = 1. - thisGal->Mgas / (thisGal->MgasIni - newMstar);
//   if(thisGal->MgasIni - newMstar == 0.)
//     feffDIVfej = 0.;
//   assert(feffDIVfej <= 1. && feffDIVfej >=0.);
  float ejectedGasFraction = 1. - thisGal->Mgas / (thisGal->MgasIni - newMstar);
  if(thisGal->MgasIni - newMstar == 0.)
    ejectedGasFraction = 0.;
  assert(ejectedGasFraction <= 1. && ejectedGasFraction >= 0.);
  if(thisGal->Mgas > 0.)
  {
    // Note that tmpMgas now is (thisGal->MgasIni - newMstar + thisGal->MgasNew)
//     MgasEj = tmpMgas * feffDIVfej;
    MgasEj = tmpMgas * ejectedGasFraction;
    MmetalEj   = MgasEj * metallicity;
    MmetalEjO  = MgasEj * metallicityO;
    MmetalEjFe = MgasEj * metallicityFe;
    MdustEj    = MgasEj * dustfraction;

    // Metal loading factor for ejection
    zetaej  = simParam->metal_ejectLoadingFactor;
    /*
    zetaej        = 2 + pow((8 / log10(thisGal->Mvir)),5);
    zetaej        = 1.0;//0.5 + pow((1e7 / thisGal->Mvir),0.3);
    */
    tmpMgas = tmpMgas - MgasEj;
    thisGal->Mmetal[0] = thisGal->Mmetal[0] - (zetaej * MmetalEj);
    thisGal->Mmetal[1] = thisGal->Mmetal[1] - (zetaej * MmetalEjO);
    thisGal->Mmetal[2] = thisGal->Mmetal[2] - (zetaej * MmetalEjFe);
    thisGal->Mdust     = thisGal->Mdust  - (zetaej * MdustEj);
  }
  else
  {
    MgasEj = tmpMgas;
    MmetalEj   = thisGal->Mmetal[0];
    MmetalEjO  = thisGal->Mmetal[1];
    MmetalEjFe = thisGal->Mmetal[2];
    MdustEj    = thisGal->Mdust;

    tmpMgas = 0.;
    thisGal->Mmetal[0] = 0.;
    thisGal->Mmetal[1] = 0.;
    thisGal->Mmetal[2] = 0.; 
    thisGal->Mdust  = 0.;
  }
  if(thisGal->Mmetal[0] < 0.) thisGal->Mmetal[0] = 0.;
  if(thisGal->Mmetal[1] < 0.) thisGal->Mmetal[1] = 0.;
  if(thisGal->Mmetal[2] < 0.) thisGal->Mmetal[2] = 0.;
  if(thisGal->Mdust < 0.) thisGal->Mdust = 0.;

  // Save metal ejected to output
  thisGal->MgasEj = MgasEj;
  thisGal->MmetalEj[0] = MmetalEj;
  thisGal->MmetalEj[1] = MmetalEjO;
  thisGal->MmetalEj[2] = MmetalEjFe;
  thisGal->MdustEj = MdustEj;
 
  thisGal->Mgas = tmpMgas;

  /* storing metal mass history of the galaxy */
  if(thisGal->Mgas != 0) 
    track_metalmasshistory(thisGal, descGal, simParam, thisGal->Mmetal[0]/thisGal->Mgas);
  else 
    track_metalmasshistory(thisGal, descGal, simParam, 0.);
  
//   if(thisGal->Mvir > 5.e9) printf("Mvir = %e\t Mdust = %e\t MdustEj = %e\n", thisGal->Mvir, thisGal->Mdust, MdustEj);
}
#endif

/*---------------------------------------------------------------------*/
/* Copy stellar, gas and DM to descendant galaxy (required to compute  */
/* accreted DM and gas masses and merged gas mass)                     */
/*---------------------------------------------------------------------*/
void copy_properties_to_descGal(dconfObj_t simParam, gal_t *thisGal, gal_t *descGal)
{
  if((descGal->snapnumber != simParam->endSnap) && (thisGal->localDescID != -1))
  {
    descGal->Mvir_prog += thisGal->Mvir;
    descGal->Mgas      += thisGal->Mgas;
    descGal->Mstar     += thisGal->Mstar;
    
#if defined WITHMETALS
    if(simParam->metals == 1)
    {
      descGal->Mmetal[0]    += thisGal->Mmetal[0];
      descGal->Mmetal[1]    += thisGal->Mmetal[1];
      descGal->Mmetal[2]    += thisGal->Mmetal[2];
      descGal->Mdust     += thisGal->Mdust;
    }
#endif
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