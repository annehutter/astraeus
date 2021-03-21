#ifndef STATISTICS_EVOLUTION_ANALYTIC_H
#define STATISTICS_EVOLUTION_ANALYTIC_H

void calc_1D_histogram_evolution_analytic(int numGal, double *property, double *numDensCurrSnap, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int minIndex1, int containsHistories1, int numBins1, int numSnaps, float **histogram, int **histogramNum);

#endif