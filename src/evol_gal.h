#ifndef EVOL_GAL_H
#define EVOL_GAL_H

void evolve_gal(gal_t *thisGal, gtree_t *thisGtree, dconfObj_t simParam);

void start_stellarmasshistory(gal_t *thisGal);
void track_stellarmasshistory(gal_t *thisGal, gal_t *descGal, dconfObj_t simParam, float newMstar);
void clean_gal(gal_t *thisGal, int outSnap);

#endif