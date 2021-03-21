#ifndef NUMDENS_ANALYTIC_H
#define NUMDENS_ANALYTIC_H

char *create_filename_hmf(char *filename, double redshift);
void read_HMF_file(char *fileName, int *numBins, double **halomass, double **hmf, double hubble_h);
double calc_numDens(int numBins, double *halomass, double *hmf, double thisHalomass, double thisBinwidth);
double *get_numDens(int numGal, double *halomassEndSnap, double thisBinwidth, char *fileName, double hubble_h);

#endif