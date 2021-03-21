#ifndef COMM_GAL_GRID_STRUCT_H
#define COMM_GAL_GRID_STRUCT_H

typedef struct{
    int32_t numGal;             /* number of galaxies allocated */
    int32_t numGalWritten;      /* number of galaxies actually written */
    float *zreion;
    float *photHI;
    int32_t *pos;
} commGalGrid_t;

commGalGrid_t *initCommGalGrid(int32_t numGal);
commGalGrid_t *initCommGalGrid_pos(int32_t numGal);
commGalGrid_t *initCommGalGrid_(int32_t numGal);
void allocCommGalGrid_zreion_photHI(commGalGrid_t **thisCommGalGrid, int32_t numGal);
void reallocCommGalGrid(commGalGrid_t **thisCommGalGrid, int32_t numGal);
void reallocCommGalGrid_pos(commGalGrid_t **thisCommGalGrid, int32_t numGal);
void reallocCommGalGrid_(commGalGrid_t **thisCommGalGrid, int32_t numGal);
void deallocate_commGalGrid(commGalGrid_t *thisCommGalGrid);

#endif