#ifndef FESC_H
#define FESC_H

double get_fesc(gal_t *thisGal, dconfObj_t simParam);
double get_fesc_const(dconfObj_t simParam);
double get_fesc_MHinc(gal_t *thisGal, dconfObj_t simParam);
double get_fesc_MHdec(gal_t *thisGal, dconfObj_t simParam);
double get_fesc_SN(gal_t *thisGal, dconfObj_t simParam);

#endif