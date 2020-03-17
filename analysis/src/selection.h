#ifndef SELECTION_H
#define SELECTION_H

void get_2D_histogram_history(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, char *filename, int thisRank);
void get_1D_histogram_history(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, int binsInLog1, int binsPerMag1, int containsHistories1, char *filename, int thisRank);

void get_2D_histogram(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, char *filename, int thisRank);
void get_1D_histogram(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *binProperty1Name, int binsInLog1, int binsPerMag1, int containsHistories1, int cumulative, char *filename, int thisRank);

char *create_filename_2D_histogram(char *directory, char *property, char *binProperty1, char *binProperty2, float redshift);
char *create_filename_1D_histogram(char *directory, char *property, char *binProperty, float redshift);
char *create_filename_2D_histogram_numDens(char *directory, char *binProperty1, char *binProperty2, float redshift);
char *create_filename_1D_histogram_numDens(char *directory, char *property, float redshift);

#endif