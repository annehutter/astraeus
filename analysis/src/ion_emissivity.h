#ifndef ION_EMISSIVITY_H
#define ION_EMISSIVITY_H

double get_nion_sps(double *stellarmasshistory, int currSnap, dconfObj_t simParam);
double get_nion_cont(double *stellarmasshistory, int currSnap, double *corrFactor_nion);
double get_nion_S99(double *stellarmasshistory, int currSnap, float *times);
double nion_S99(float age);
double get_nion_BPASS(double *stellarmasshistory, int currSnap, float *times);
double nion_BPASS(float age);

double *get_corrFactor_nion(dconfObj_t simParam, int endSnap);
double *get_corrFactor_nion_S99(int endSnap, float *times);
double *get_corrFactor_nion_BPASS(int endSnap, float *times);
double *create_corrFactor_nion(int endSnap, float *times, float t0, float tbreak, float exponent, double nion);
float calc_corrFactor_nion(float tinitial, float tfinal, float t0, float tbreak, float exponent);

#endif