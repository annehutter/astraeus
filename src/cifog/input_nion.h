#ifndef INPUT_NION_H
#define INPUT_NION_H

void read_update_nion(confObj_t simParam, fftw_complex *nion, grid_t *thisGrid);
void read_update_nion_HeI(confObj_t simParam, fftw_complex *nion_HeI, grid_t *thisGrid);
void read_update_nion_HeII(confObj_t simParam, fftw_complex *nion_HeII, grid_t *thisGrid);

#endif