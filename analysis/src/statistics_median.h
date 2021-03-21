#ifndef STATISTICS_MEDIAN_H
#define STATISTICS_MEDIAN_H

void calc_2D_histogram_history_median(int numGal, double *property, int numSnaps, float *redshifts, int currSnap, double *binProperty1, double *binProperty2, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, int size, char *filename);
void calc_2D_histogram_median(int numGal, double *property, int currSnap, double *binProperty1, double *binProperty2, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, int size, char *filename);

void write_2D_histogram_history_median(int numBins1, float *values1, int numBins2, float *values2, int currSnap, float *redshifts, float *histogram, int *numHistogram, int thisRank, char *filename);
void write_2D_histogram_median(int numBins1, float *values1, int numBins2, float *values2, float *histogram, int *numHistogram, int thisRank, char *filename);

#endif