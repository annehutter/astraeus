#ifndef READ_SNAP_H
#define READ_SNAP_H

outgalsnap_t *read_snap(int snap, int nb_gal, char *baseName);
outgalsnap_t **read_snaps(int minSnap, int maxSnap, int *nb_gals);
int *get_nblines();

#endif
