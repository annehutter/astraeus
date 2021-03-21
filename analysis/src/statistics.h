#ifndef STATISTICS_H
#define STATISTICS_H

int get_max(int *arrayInt, int length);
void get_min_and_max_galaxy_property(int numGal, double *property, float *min, float *max);
void get_min_and_max_galaxy_property_with_histories(int numGal, double *property, float *min, float *max, int currSnap);
void get_bins(float minValue, float maxValue, int binsInLog, int binsPerMag, int *numBinsIn, float **valuesIn, int *minIndexIn);

void calc_1D_histogram_history(int numGal, double *property, int numSnaps, float *redshifts, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1, int thisRank, char *filename);
void calc_2D_histogram_history(int numGal, double *property, int numSnaps, float *redshifts, int currSnap, double *binProperty1, double *binProperty2, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, char *filename);

void write_1D_histogram_history(int numBins1, float *values1, int currSnap, float *redshifts, float *histogram, int *numHistogram, int thisRank, char *filename);
void write_2D_histogram_history(int numBins1, float *values1, int numBins2, float *values2, int currSnap, float *redshifts, float *histogram, int *numHistogram, int thisRank, char *filename);

void calc_1D_histogram_value(int numGal, double *property, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1, int thisRank, char *filename);
void calc_2D_histogram_value(int numGal, double *property, int currSnap, double *binProperty1, double *binProperty2, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, char *filename);
void calc_3D_histogram_value(int numGal, double *property, int currSnap, double *binProperty1, double *binProperty2, double *binProperty3, int binsInLog1, int binsInLog2, int binsInLog3, int binsPerMag1, int binsPerMag2, int binsPerMag3, int containsHistories1, int containsHistories2, int containsHistories3, int thisRank, char *filename);

void write_1D_histogram_value(int numBins1, float *values1, float *histogram, int *numHistogram, int thisRank, char *filename);
void write_2D_histogram_value(int numBins1, float *values1, int numBins2, float *values2, float *histogram, int *numHistogram, int thisRank, char *filename);
void write_3D_histogram_value(int numBins1, float *values1, int numBins2, float *values2, int numBins3, float *values3, float *histogram, int *numHistogram, int thisRank, char *filename);

void calc_1D_histogram(int numGal, double *binProperty1, int currSnap, int binsInLog1, int binsPerMag1, int containsHistories1, float volumeInMpc, int cumulative, int thisRank, char *filename);
void calc_2D_histogram(int numGal, double *binProperty1, double *binProperty2, int currSnap, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, float volumeInMpc, int thisRank, char *filename);

void write_1D_histogram(int num_bins, float *values, double *hist, int thisRank, char *filename);
void write_2D_histogram(int numBins1, float *values1, int numBins2, float *values2, double *numHistogram, int thisRank, char *filename);

#endif
