#ifndef METALLICITY_H
#define METALLICITY_H

#if defined WITHMETALS

typedef struct{
  int32_t endSnap;
  double *meanMetallicityIgm;
  double *meanMetallicityIgmO;
  double *meanMetallicityIgmFe;
  double *meanDustfractionIgm;
} metal_t;

metal_t *initMetal(int32_t endSnap);
void deallocate_metal(metal_t *thisMetal);

int read_array_from_txtfile(char *fileName, double **thisArray);

int where_loop(double array[], int n, double value);
int where(double array[], int n, double search);
double linear_interp(double x[], double y[], int length, double xx);

void add_ejected_metals_to_igm(gal_t *thisGal, int snap, double invTotMgas, metal_t *thisMetal);
void add_igm_metals_to_gal(metal_t *thisMetal, gal_t *thisGal, int snap);
void map_metallicityIgm_to_gal(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, metal_t *thisMetal);
void distribute_metallicityIgm(metal_t *thisMetal, int snap);

int *compute_indexes(int lower, int upper, int numMetallicityBins, double *metallicity, double *metallicityAtFormation);
double compute_yield(double *yields, double *numStarsAtFormation, int numMassBins, int *indexes, int index_low, int index_up, int lowerIntegrationLimit, int upperIntegrationLimit, double stepsizeIntegration);
double compute(double yield, double *numStarsAtFormation, int lowerIntegrationLimit, int upperIntegrationLimit, double stepsizeIntegration);

void compute_metal_mass(gal_t *thisGal, dconfObj_t simParam, float sfr, float deltat, float *G_t, float *e_Z, float *e_Z_O, float *e_Z_Fe, float *y_d);
#endif 

#endif
