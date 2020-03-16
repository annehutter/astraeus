#ifndef ION_EMISSIVITY_H
#define ION_EMISSIVITY_H

double get_nion_for_model(gal_t *thisGal, dconfObj_t simParam);
double get_nion_sps(gal_t *thisGal, dconfObj_t simParam);
double get_nion_cont(gal_t *thisGal, double *corrFactor_nion);
double get_nion_S99(gal_t *thisGal, float *times);
double nion_S99(float age);
double get_nion_BPASS(gal_t *thisGal, float *times);
double nion_BPASS(float age);

double *get_corrFactor_nion(dconfObj_t simParam);
double *get_corrFactor_nion_S99(int endSnap, float *times);
double *get_corrFactor_nion_BPASS(int endSnap, float *times);
double *create_corrFactor_nion(int endSnap, float *times, float t0, float tbreak, float exponent, double nion);
float calc_corrFactor_nion(float tinitial, float tfinal, float t0, float tbreak, float exponent);

#endif