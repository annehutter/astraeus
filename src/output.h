#ifndef OUTPUT_H
#define OUTPUT_H

/*-------------------------------------------------------*/
/* COPYING TO OUTPUT STRUCTS                             */
/*-------------------------------------------------------*/

/* VERTICAL OUPUT */
void copy_tree_to_outgtree(gtree_t *thisTree, outgtree_t *thisOutGTree);
void copy_gal_to_outgaltree(gal_t *thisGal, outgal_t *thisOutGal);

/* HORIZONTAL OUPUT */
void copy_gal_to_outgalsnap(gal_t *thisGal, outgalsnap_t *thisOutGal);
void copy_gal_to_outGalList(gal_t *thisGal, int thisOutGal, outgalsnap_t *thisOutGalList);

/*-------------------------------------------------------*/
/* WRITING OUTPUT FILES                                  */
/*-------------------------------------------------------*/

char *define_outfilename(int snap, const char *genericOutfileName );

/* VERTICAL OUPUT */
void write_treelist(dconfObj_t simParam, int thisRank, int numGTree, gtree_t **thisGtreeList);

/* HORIZONTAL OUPUT */
void write_galaxies_of_snap_to_file(dconfObj_t simParam, int snap, int numOutGal, outgalsnap_t *outGalList);
void write_num_galaxies_to_file(dconfObj_t simParam, int *numGalSnap, int numSnaps, int thisRank);

/*-------------------------------------------------------*/
/* DELETING EXISTING FILES                               */
/*-------------------------------------------------------*/

void delete_file(int snap, const char *genericOutfileName);
void delete_allfiles(int startSnap, int endSnap, int deltaSnap, const char *genericOutfileName);

/*-------------------------------------------------------*/
/* COPYINH INIFILE                                       */
/*-------------------------------------------------------*/

void copy_iniFile_executable(char *outputfile);

#endif
