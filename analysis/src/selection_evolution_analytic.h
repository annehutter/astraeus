#ifndef SELECTION_EVOLUTION_ANALYTIC_H
#define SELECTION_EVOLUTION_ANALYTIC_H

void get_1D_histogram_evolution_analytic(dconfObj_t simParam, outgtree_t **theseTrees, int numTrees, int **index, int ***listEquals, int *sizeListEquals, int **listMerged, int listID, int thisRank);

#endif