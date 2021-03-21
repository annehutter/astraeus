#ifndef OUTPUT_LISTS_H
#define OUTPUT_LISTS_H

void write_galaxies_to_textfile_all_snaps(dconfObj_t simParam, outgtree_t **theseTrees, int numTrees, int **index, int ***listEquals, int *sizeListEquals, int **listMerged, int thisRank);
void write_galaxies_to_textfile(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *filename, int thisRank);
char *create_filename_output_lists(char *directory, float redshift, int thisRank);

#endif