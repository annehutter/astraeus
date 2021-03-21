#ifndef PHOTION_GAL_H
#define PHOTION_GAL_H

double get_gridProperty_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *gridProperty);
void update_field(fftw_complex **field, grid_t *thisGrid, double redshift, double XHII_threshold, char *fieldName, int memoryIntensive);
void update_grid_field(grid_t *thisGrid, double redshift, double XHII_threshold, char *fieldName);
void update_local_field(fftw_complex **field, grid_t *thisGrid, double XHII_threshold, char *fieldName);

void write_grid_to_file_double_test(fftw_complex *thisArray, int nbins, int local_n0, int local_0_start, char *filename, int thisRank, int memory_intensive);

#endif
