#ifndef RUN_H
#define RUN_H

void run_astraeus(dconfObj_t simParam, int thisRank, int size);
#if defined WITHMETALS
void evolve(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, dconfObj_t simParam, domain_t *thisDomain, nion_t *thisNionList, fftw_complex *photHI, fftw_complex *zreion, metal_t *thisMetal, int outSnap, outgalsnap_t **outGalList, int *numGalInSnapList);
#else
void evolve(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, dconfObj_t simParam, domain_t *thisDomain, nion_t *thisNionList, fftw_complex *photHI, fftw_complex *zreion, int outSnap, outgalsnap_t **outGalList, int *numGalInSnapList);
#endif
void set_walker_to_startSnap(int numGtrees, gtree_t ***thisGtreeList, int startSnap);

#if defined WITHMETALS
void get_and_initialise_metallicity_tables(dconfObj_t simParam);
#endif

void write_zreion_grids_to_file(dconfObj_t simParam, domain_t *thisDomain, int snap, grid_t *thisGrid, fftw_complex *zreion, fftw_complex *photHI);

#endif

