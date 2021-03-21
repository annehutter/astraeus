#ifndef SELECTION_EVOLUTION_H
#define SELECTION_EVOLUTION_H

void get_1D_histogram_evolution(dconfObj_t simParam, outgtree_t **theseTrees, int numTrees, int **index, int ***listEquals, int *sizeListEquals, int **listMerged, int listID, int thisRank);
char *create_filename_1D_histogram_evolution(char *directory, char *property, char *binProperty, float smoothingScale, double MvirThreshold);

#endif