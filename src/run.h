#ifndef RUN_H
#define RUN_H

void run_delphi(dconfObj_t simParam, int thisRank, int size);
void evolve(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, dconfObj_t simParam, domain_t *thisDomain, nion_t *thisNionList, fftw_complex *photHI, fftw_complex *XHII, fftw_complex *zreion, int outSnap, outgalsnap_t **outGalList, int *numGalInSnapList);
void set_walker_to_startSnap(int numGtrees, gtree_t ***thisGtreeList, int startSnap);

#endif

