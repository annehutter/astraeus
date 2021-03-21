#ifndef STATISTIC_ANALYTIC_H
#define STATISTIC_ANALYTIC_H

void calc_1D_histogram_analytic(int numGal, double *binProperty1, double *numDensEndSnap, int currSnap, int binsInLog1, int binsPerMag1, int containsHistories1, int cumulative, int thisRank, char *filename);

void calc_2D_histogram_analytic(int numGal, double *binProperty1, double *binProperty2, double *numDensEndSnap, int currSnap, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, char *filename);

void write_1D_histogram_analytic(int num_bins, float *values, double *hist, int thisRank, char *filename);

void write_2D_histogram_analytic(int numBins1, float *values1, int numBins2, float *values2, double *numHistogram, int thisRank, char *filename);

#endif