#ifndef GET_GRID_PROPERTIES_H
#define GET_GRID_PROPERTIES_H

/*----------------------------------------------------------------*/
/* DO DOMAIN DECOMPOSITION FOR GRIDS */
/*----------------------------------------------------------------*/

// int get_grid_decomposition(int size, int thisRank, int gridsize);

/*----------------------------------------------------------------*/
/* READING GRIDS */
/*----------------------------------------------------------------*/

float *read_grid(char *filename, int gridsize, int doublePrecision, int memoryIntensive);
void read_grid_singleprecision(float *toThisArray, int gridsize, int zNbins, int thisOffset, char *filename);
void read_grid_doubleprecision(float *toThisArray, int gridsize, int zNbins, int thisOffset, char *filename);

/*----------------------------------------------------------------*/
/* COMPUTE XHI GRID FROM XHII GRID */
/*----------------------------------------------------------------*/

void get_XHI_grid(float **thisGrid, int gridsize, int memoryIntensive);

/*----------------------------------------------------------------*/
/* GETTING GRID PROPERTIES FOR GALAXIES */
/*----------------------------------------------------------------*/

int check_grid_property_present(char *propertyName);
void get_grid_values(int numEntries, double **toThisArray, double *posIndex, int gridsize, float *thisGrid, int memoryIntensive);
double *get_grid_property(int size, int thisRank, int gridProperty, char *gridFilename, int numGal, double *posIndex, int gridsize, int doublePrecision, int memoryIntensive, int smoothedGrid, double smoothingScale);
double *getGridProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times);

#endif
