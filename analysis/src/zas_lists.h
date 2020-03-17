#ifndef ZAS_LISTS_H
#define ZAS_LISTS_H

int get_num_snapshots(outgtree_t **thisTreeList);
float *get_scalefactors(outgtree_t **thisTreeList, int numTrees);
float *get_redshifts(float *scalefactors, int numSnaps);
float *get_times_from_redshifts(float *redshifts, int numSnaps, double h, double omega_m, double omega_l);
float calc_time_from_redshift(double zmin, double zmax, double h, double omega_m, double omega_l);
int find_snapshoft_from_redshift(float redshift, int numSnaps, float *redshifts);

#endif