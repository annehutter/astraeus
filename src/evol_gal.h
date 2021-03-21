#ifndef EVOL_GAL_H
#define EVOL_GAL_H

void evolve_gal(gal_t *thisGal, gtree_t *thisGtree, dconfObj_t simParam);

void start_stellarmasshistory(gal_t *thisGal);
void track_stellarmasshistory(gal_t *thisGal, gal_t *descGal, dconfObj_t simParam, float newMstar);
void clean_gal(gal_t *thisGal, int outSnap);

float calc_Rvir(float Mvir, double omega_m);
void correct_Mvir_to_continously_rise(gal_t *thisGal);
void get_accreted_and_merged_gas(dconfObj_t simParam, gal_t *thisGal);
void apply_radFeedback_and_update_zreion(dconfObj_t simParam, gal_t *thisGal, gal_t * descGal);
void store_available_gas(gal_t *thisGal);
void do_starFormation_and_SNfeedback(dconfObj_t simParam, gal_t *thisGal, gal_t *descGal);
void copy_properties_to_descGal(dconfObj_t simParam, gal_t *thisGal, gal_t *descGal);

void pass_zreion_to_descGal(gal_t *thisGal, gal_t *descGal, dconfObj_t simParam);

#endif