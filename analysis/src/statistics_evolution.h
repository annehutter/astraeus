#ifndef STATISTICS_EVOLUTION_H
#define STATISTICS_EVOLUTION_H

int get_numBins_for_bins(int numGal, int currSnap, float *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1);
int get_minIndex_for_bins(int numGal, int currSnap, float *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1);
float *get_values_for_bins(int numGal, int currSnap, float *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1);
int get_numBins_for_fixed_bins(float minValue1, float maxValue1, int binsInLog1, int binsPerMag1);
int get_minIndex_for_fixed_bins(float minValue1, float maxValue1, int binsInLog1, int binsPerMag1);
float *get_values_for_fixed_bins(float minValue1, float maxValue1, int binsInLog1, int binsPerMag1);

void calc_1D_histogram_evolution(int numGal, float *property, int currSnap, float *binProperty1, int binsInLog1, int binsPerMag1, int minIndex1, int containsHistories1, int numBins1, int numSnaps, float **histogram, int **histogramNum);

void write_1D_histogram_evolution(int numBins, float *values, int endSnap, float *redshifts, float *histogram, int *numHistogram, int thisRank, char *filename);

#endif