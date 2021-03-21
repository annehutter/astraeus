#ifndef SELECTION_MEDIAN_H
#define SELECTION_MEDIAN_H

void get_2D_histogram_history_median(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, char *filename, int thisRank, int size);
void get_2D_histogram_median(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *propertyName, char *binProperty1Name, char *binProperty2Name, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, double binProperty1LowLimit, double binProperty1UpLimit, char *filename, int thisRank, int size);

char *create_filename_2D_histogram_history_median(char *directory, char *property, char *binProperty1, char *binProperty2, float redshift, float smoothingScale, double MvirThreshold);
char *create_filename_2D_histogram_median(char *directory, char *property, char *binProperty1, char *binProperty2, float redshift, float smoothingScale, double MvirThreshold);

#endif