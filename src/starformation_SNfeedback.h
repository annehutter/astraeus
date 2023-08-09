#ifndef STARFORMATION_SNFEEDBACK_H
#define STARFORMATION_SNFEEDBACK_H

/*---------------------------------------------------------------------*/
/* SUPERNOVA FEEDBACK                                                  */
/*---------------------------------------------------------------------*/

float calc_SNejection_efficiency(gal_t *thisGal, dconfObj_t simParam, float corrFactorTimeStep);
float calc_remaining_gasfraction(gal_t *thisGal, dconfObj_t simParam, float corrFactorTimeStep);
float get_SNenergy_current(gal_t *thisGal, float *SNenergy);
float calc_SNenergy_past(gal_t *thisGal, float *SNenergy);

/*---------------------------------------------------------------------*/
/* COMPUTING DELAYED SUPERNOVA FEEDBACK TABLE                          */
/*---------------------------------------------------------------------*/

float *get_SNenergy_delayed(dconfObj_t simParam);
float get_SNfraction_per_Msun(float tsnap, float tsnapPrev, float tprevSnap, float tprevprevSnap, int burstySF);
float get_SNmass(float timeInMyr);
float get_SNtime(float massInMsun);
float calc_SNfraction_per_Msun(float MSNlow, float MSNup, float Mlow, float Mup, float slope);
float calc_SNfraction_per_Msun_contSF(float tsnap, float tsnapPrev, float tprevSnap, float tprevprevSnap, float MSNlow, float Mlow, float Mup, float slope);

#endif