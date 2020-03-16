#ifndef RADIATIVE_FEEDBACK_H
#define RADIATIVE_FEEDBACK_H

float calc_radiative_feedback(gal_t *thisGal, gal_t *descGal, dconfObj_t simParam);

/*---------------------------------------------------------------------*/
/* MASS LIMITS                                                         */
/*---------------------------------------------------------------------*/
float get_Mcool(float Temp, float mu, float z, float omega_m);
float get_Mjeans(float Temp, float mu, float z, float omega_m);
float get_Mjeans_vir(float Temp, float mu, float z, float omega_m);

/*---------------------------------------------------------------------*/
/* EXTREME RADIATIVE FEEDBACK DESCRIPTION                              */
/*---------------------------------------------------------------------*/
float fg_radfb_Mjeans(float Mvir, float Temp, float mu, float z, float omega_m);
float fg_radfb_Mcool(float Mvir, float Temp, float mu, float z, float omega_m);

/*---------------------------------------------------------------------*/
/* RADIATIVE FEEDBACK DESCRIPTION FOLLOWING SOBACCHI 2013              */
/*---------------------------------------------------------------------*/
float fg_radfb_Sobacchi(float photHI, float Mvir, float z, float zreion);

/*---------------------------------------------------------------------*/
/* RADIATIVE FEEDBACK DESCRIPTION FOLLOWING GNEDIN 2000                */
/*---------------------------------------------------------------------*/
float fg_radfb_Gnedin(float Mvir, float Temp, float mu, float z, float zion, float zreion, float omega_m);
float fg_radfb_TempEvol(float Mvir, float Temp, float mu, float z, float zion, float zreion, float omega_m);

float get_fa_Gnedin(float z, float zion, float zreion);
float get_fa(float z, float zion, float zreion, float temp);

float get_frec(float a, float ap, float arec, float temp);
float get_fdec(float a, float ap, float adec, float temp);
float get_fheat(float a, float ap, float aion, float aheat, float alpha);
float get_fion(float a, float ap, float aion);
float get_fcool(float a, float areion);

#endif
