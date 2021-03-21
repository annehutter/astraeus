#ifndef NION_GAL_H
#define NION_GAL_H

/*-----------------------------------------------------*/
/* MAPPING GALAXIES ONTO GRID                          */
/*-----------------------------------------------------*/

/* GENERATE FFTW ARRAY FOR NION */
fftw_complex *map_galnion_to_grid(nion_t *thisNion, domain_t *thisDomain, int memoryIntensive);

/* MAPPING FUNCTION TO GRID */
void map_nion_to_pos_on_grid(nion_t *thisNion, fftw_complex *nion, int nbins, ptrdiff_t local_n0, ptrdiff_t local_0_start, int thisRank);

/* MEMORY INTENSIVE MAPPING FUNCTION TO GRID */
void map_nion_to_pos_on_grid_memintensive(nion_t *thisNion, fftw_complex *nion, int nbins, ptrdiff_t local_n0, ptrdiff_t local_0_start, int thisRank);

/*=======================================================================================*/

/*-----------------------------------------------------*/
/* COMMUNICATION FOR NION OF GALAXIES TO BE MAPPED     */
/*-----------------------------------------------------*/

/* MAIN ROUTINE: SEND & RECIVE GALAXIES TO PROCESSORS ACCORDING TO THEIR LOCATION */
#ifdef MPI
void distribute_nion_to_processors(nion_t *thisNion, domain_t *thisDomain);
#endif

/* SORT ARRAY ACCORDING TO RANK WHERE TO SEND IT */
void sort_nion(nion_t *thisNion, int thisRank);

/* GET HOW MANY GALAXIES HAVE TO BE SENT TO EACH RANK (works only on a sorted array!) */
int32_t *get_nion_numbers(nion_t *thisNion, int32_t size);

/* GENERATE A COPY OF NION IN THISNION */
double *copy_nion(nion_t *thisNion);

/* GENERATE A COPY OF POS IN THISNION */
int32_t *copy_pos(nion_t *thisNion);

/* SEND & RECEIVE NION & POS TO PROCESSORS ACCORDING TO THEIR LOCATION */
#ifdef MPI
void send_recv_nion_pos(domain_t *thisDomain, int32_t *numGalToRanks, int32_t **numGalOnThisRank, double *toSendNion, int32_t *toSendPos, double **toRecvNion, int32_t **toRecvPos);
#endif

/*=======================================================================================*/

/*-----------------------------------------------------*/
/* COPYING GALAXY PROPERTIES TO NION LIST TO BE MAPPED */
/*-----------------------------------------------------*/

void copy_gal_to_nion(gal_t *thisGal, dconfObj_t simParam, domain_t *thisDomain, double factor, nion_t *thisNion);

#endif