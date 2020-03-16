#ifndef REDSHIFT_LIST_H
#define REDSHIFT_LIST_H

typedef struct{
  int numRedshifts;
  double *redshifts;
}redshiftlist_t;

int get_num_lines(char *fileName);
redshiftlist_t *read_redshift_snap_list(char *fileName);
void deallocate_redshiftlist(redshiftlist_t *thisRedshiftList);
double calc_deltaRedshift(redshiftlist_t *thisRedshiftList, int snap);
int get_snap_from_redshift(redshiftlist_t *thisRedshiftList, double redshift);

float *get_times_from_redshifts(redshiftlist_t *thisRedshiftList, int num, double h, double omega_m, double omega_l);
float calc_time_from_redshift(double zmin, double zmax, double h, double omega_m, double omega_l);

#endif