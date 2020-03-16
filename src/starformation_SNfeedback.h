#ifndef STARFORMATION_SNFEEDBACK_H
#define STARFORMATION_SNFEEDBACK_H

/*---------------------------------------------------------------------*/
/* SUPERNOVA FEEDBACK                                                  */
/*---------------------------------------------------------------------*/

float calc_SNejection_efficiency(gal_t *thisGal, dconfObj_t simParam, float Mgas);
float get_SNenergy_current(gal_t *thisGal, float *SNenergy);
float calc_SNenergy_past(gal_t *thisGal, float *SNenergy);

/*---------------------------------------------------------------------*/
/* COMPUTING DELAYED SUPERNOVA FEEDBACK TABLE                          */
/*---------------------------------------------------------------------*/

float *get_SNenergy_delayed(int endSnap, float *times);
float get_SNfraction_per_Msun(float tsnap, float tprevSnap, float tprevprevSnap);
float get_SNmass(float timeInGyr);
float calc_SNfraction_per_Msun(float MSNlow, float MSNup, float Mlow, float Mup, float slope);

#endif