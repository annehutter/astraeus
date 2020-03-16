#ifndef CIFOG_DELPHI_H
#define CIFOG_DELPHI_H

int cifog_init(char *iniFile, confObj_t *simParam, grid_t **grid,  integral_table_t **integralTable, photIonlist_t **photIonBgList, const int myRank);

int cifog(confObj_t simParam, int snap, double redshift, double deltaRedshift, grid_t *grid, fftw_complex *nion, fftw_complex *nion_HeI, fftw_complex *nion_HeII, const integral_table_t *integralTable, photIonlist_t *photIonBgList, const int myRank);

int cifog_step(confObj_t simParam, grid_t *grid, fftw_complex *nion, fftw_complex *nion_HeI, fftw_complex *nion_HeII, const integral_table_t *integralTable, photIonlist_t *photIonBgList, const int cycle, const int cycle_offset, int snap, const int myRank);

int cifog_deallocate(confObj_t simParam, grid_t *grid, integral_table_t *integralTable, photIonlist_t *photIonBgList, const int myRank);

#endif