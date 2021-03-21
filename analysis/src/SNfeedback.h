#ifndef STARFORMATION_SNFEEDBACK_H
#define STARFORMATION_SNFEEDBACK_H

/*---------------------------------------------------------------------*/
/* SUPERNOVA FEEDBACK                                                  */
/*---------------------------------------------------------------------*/

double *get_SNejection_efficiency(dconfObj_t simParam, int numGal, double *MvirDivRvir, double *MgasIni, double *stellarMassHistory, int currSnap);
double calc_SNejection_efficiency(dconfObj_t simParam, double MvirDivRvir, double Mgas, double *stellarmasshistory, int currSnap);
float get_SNenergy_current(int currSnap, float *SNenergy);
float calc_SNenergy_past(int currSnap, float *SNenergy, double *stellarmasshistory);

/*---------------------------------------------------------------------*/
/* COMPUTING DELAYED SUPERNOVA FEEDBACK TABLE                          */
/*---------------------------------------------------------------------*/

float *get_SNenergy_delayed(int endSnap, float *times);
float get_SNfraction_per_Msun(float tsnap, float tprevSnap, float tprevprevSnap);
float get_SNmass(float timeInGyr);
float calc_SNfraction_per_Msun(float MSNlow, float MSNup, float Mlow, float Mup, float slope);

#endif