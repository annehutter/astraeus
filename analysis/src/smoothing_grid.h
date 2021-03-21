#ifndef SMOOTHING_GRID_H
#define SMOOTHING_GRID_H

void smooth_grid_property(float **thisGrid, int nbins, int memoryIntensive, float smooth_scale);
void get_grid_decomposition(int nbins, int *local_n0, int *local_0_start, int memoryIntensive);
fftw_complex *construct_tophat_filter(int nbins, ptrdiff_t local_0_start, ptrdiff_t local_n0, float smooth_scale, int memoryIntensive);
void convolve_fft(int nbins, ptrdiff_t local_0_start, ptrdiff_t local_n0, fftw_complex *filter, fftw_complex **output, fftw_complex **input, int memoryIntensive);

#endif