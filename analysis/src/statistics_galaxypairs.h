#ifndef STATISTICS_GALAXYPAIRS_H
#define STATISTICS_GALAXYPAIRS_H

void send_recv_indeces_pos_primaryGalaxies(int *numPrimaryGalOnRanks, int **primaryGalIndex, double **primaryGalPosx, double **primaryGalPosy, double **primaryGalPosz, int thisRank, int size);

double get_distanceSqared_galaxy_pair(double posx1, double posy1, double posz1, double posx2, double posy2, double posz2);

void search_secondaryGalaxies(int numPrimaryGal, int *primaryGalOriginRanks, double *primaryGalPosx, double *primaryGalPosy, double *primaryGalPosz, int numGal, double *selectionProperty, double minSelectionProperty, double maxSelectionProperty, double maxDistance, double *galPosx, double *galPosy, double *galPosz, int *galIndex, int ***galIndexAroundPrimaryGal, double ***galDistanceFromPrimaryGal, int ***galOriginRankOfPrimaryGal);

void write_galaxypairs_to_file(int numPrimaryGal, int *primaryGalIndex, int thisRank, int size, int *numPrimaryGalOnRanks, int numGal, double *Mvir, double *MgasIni, double *Mgas, double *Mstar, double *density, double *XHII, int **galIndexAroundPrimaryGal, double **galDistanceFromPrimaryGal, int **galOriginRankOfPrimaryGal, char *basename);

void write_galaxypairs_SFHs_to_file(int numPrimaryGal, int *primaryGalIndex, int thisRank, int size, int *numPrimaryGalOnRanks, int numGal, int currSnap, float *redshifts, double *propertyWithHistory, int **galIndexAroundPrimaryGal, double **galDistanceFromPrimaryGal, int **galOriginRankOfPrimaryGal, char *basename);

void statistics_on_galaxypairs(int numGal, double *selectionProperty, double minSelectionProperty, double maxSelectionProperty, double *selectionProperty2, double minSelectionProperty2, double maxSelectionProperty2, double maxDistance, int numSnaps, float *redshifts, int currSnap, double *posxProperty, double *posyProperty, double *poszProperty, double *Mvir, double *MgasIni, double *Mgas, double *Mstar, double *density, double *XHII, double *propertyWithHistory, int thisRank, int size, char *basename);

#endif