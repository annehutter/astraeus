#ifndef MAP_GAL_TO_GRID_H
#define MAP_GAL_TO_GRID_H

int check_map_gal_to_grid_property_present(char *propertyName);
char *get_map_gal_to_grid_propertyName(char *propertyName);
double *get_mapped_gal_positions(int numGal, double *galProperty, double *posIndex, double lowLimit, double upLimit, int *selectedNumGal);

#ifdef MPI
float *create_grid_with_map_gal_mpi(int size, int thisRank, int numEntries, double *posIndex, int gridsize);
float *create_grid_with_map_gal_mpi_memintensive(int numEntries, double *posIndex, int gridsize);
#endif
float *create_grid_with_map_gal_serial(int numEntries, double *posIndex, int gridsize);

float *create_grid_with_map_gal(int size, int thisRank, int numEntries, double *posIndex, int gridsize, int memoryIntensive);

double *get_map_gal_to_grid_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, double galPropertyLowLimit, double galPropertyUpLimit, int currSnap, float *times);

#endif