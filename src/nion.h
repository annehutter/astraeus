#ifndef NION_H
#define NION_H

typedef struct{
    int32_t numGal;             /* number of galaxies allocated */
    int32_t numGalWritten;      /* number of galaxies actually written */
    int32_t *rank;
    double *nion;
    int32_t *pos;
} nion_t;

nion_t *initNion(int32_t numGal);
void reallocNion(nion_t **thisNion, int32_t numGal);
void deallocate_nion(nion_t *thisNion);
double get_mean_nion(nion_t *thisNion);

#endif