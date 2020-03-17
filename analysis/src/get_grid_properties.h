#ifndef GET_GRID_PROPERTIES_H
#define GET_GRID_PROPERTIES_H

/*----------------------------------------------------------------*/
/* DO DOMAIN DECOMPOSITION FOR GRIDS */
/*----------------------------------------------------------------*/

int get_grid_decomposition(int size, int thisRank, int gridsize);

/*----------------------------------------------------------------*/
/* READING GRIDS */
/*----------------------------------------------------------------*/

float *read_grid(int size, int thisRank, char *filename, int gridsize, int doublePrecision, int memoryIntensive);
void read_grid_singleprecision(float *toThisArray, int gridsize, int zNbins, int thisOffset, char *filename);
void read_grid_doubleprecision(float *toThisArray, int gridsize, int zNbins, int thisOffset, char *filename);

/*----------------------------------------------------------------*/
/* GETTING GRID PROPERTIES FOR GALAXIES */
/*----------------------------------------------------------------*/

int check_grid_property_present(char *propertyName);
float *get_grid_property(int size, int thisRank, char *gridFilename, int numGal, float *posIndex, int gridsize, int doublePrecision, int memoryIntensive);
void get_grid_values(int size, int thisRank, int numEntries, float **toThisArray, float *posIndex, int gridsize, float *thisGrid);

#endif
