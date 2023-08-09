#ifndef SELECTION_GALAXYPAIRS_H
#define SELECTION_GALAXYPAIRS_H

void select_galaxypairs(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *selectionPropertyName, char *selectionProperty2Name, char *propertyWithHistoryName, double minSelectionProperty, double maxSelectionProperty, double minSelectionProperty2, double maxSelectionProperty2, double maxDistance, char *filename, int thisRank, int size);

char *create_filename_galaxypairs(char *directory, char *selectionPropertyName, char *selectionProperty2Name, double minSelectionProperty, double maxSelectionProperty, double minSelectionProperty2, double maxSelectionProperty2, float redshift);

#endif