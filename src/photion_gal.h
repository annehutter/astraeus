#ifndef PHOTION_GAL_H
#define PHOTION_GAL_H

void update_photHI_field(fftw_complex **photHI, grid_t *thisGrid);
void update_photHI_zreion_field(fftw_complex **photHI_zreion, grid_t *thisGrid, double redshift, double XHII_threshold);
double get_photHI_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *photHI);

void update_XHII_field(fftw_complex **XHII, grid_t *thisGrid);
double get_XHII_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *XHII);

void update_zreion_field(fftw_complex **zreion, grid_t *thisGrid, double redshift, double XHII_threshold);
double get_zreionField_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *zreion);
double get_zreion_gal(gal_t *thisGal, int32_t nbins, double inv_boxsize, fftw_complex *XHII, double XHII_threshold);

void write_grid_to_file_double_test(fftw_complex *thisArray, int nbins, int local_n0, char *filename);

#endif
