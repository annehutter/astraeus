#ifndef COMM_GAL_GRID_H
#define COMM_GAL_GRID_H

/* ------------------------------------------------------------------------*/
/* MAIN FUNCTION TO GET ZREION & PHOTHI FOR GALAXIES IN NEXT SNAPSHOT      */
/* ------------------------------------------------------------------------*/

void map_grid_to_gal(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, dconfObj_t simParam, domain_t *thisDomain, grid_t *thisGrid);

/* ------------------------------------------------------------------------*/
/* FUNCTIONS TO GET POSITIONS AND GRID RANKS FOR GALAXIES IN NEXT SNAPSHOT */
/* ------------------------------------------------------------------------*/

int32_t get_gal_position_index(gal_t *thisGal, dconfObj_t simParam, domain_t *thisDomain);
int32_t get_gal_gridRank(gal_t *thisGal, dconfObj_t simParam, domain_t *thisDomain);
void get_galaxies_next_snapshot(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, dconfObj_t simParam, domain_t *thisDomain, commGalGrid_t *thisCommGalGrid, int32_t **gridRanks);

/* ------------------------------------------------------------------------*/
/* FUNCTIONS TO GET GRID PROPERTIES FOR GALAXIES IN NEXT SNAPSHOT */
/* ------------------------------------------------------------------------*/

void get_grid_properties_galaxies_next_snapshot(domain_t *thisDomain, grid_t *thisGrid, commGalGrid_t *thisCommGalGrid, int32_t **gridRanks);
int32_t *create_index_array(int numGal);
void sort_gal_positions(commGalGrid_t *thisCommGalGrid, int32_t *indexArray);
int32_t *get_numGal_gridRanks(int numGal, int32_t *gridRanks, int32_t size);

void send_recv_pos(domain_t *thisDomain, int32_t *numGal_gridRanks, int32_t **numGal_galRanks, int32_t *toSendPos, int32_t **toRecvPos);

int32_t sum_numGal_ranks(int size, int32_t *numGal_ranks);

void send_recv_zreion_photHI(domain_t *thisDomain, int32_t *numGal_gridRanks, int32_t **numGal_galRanks, float *toSendZreion, float *toSendPhotHI, float **toRecvZreion, float **toRecvPhotHI);

int32_t check_numGal_gridRanks(int size, int32_t *numGal_gridRanks, int32_t *recvNumGal_gridRanks);
void get_gal_zreion_photHI(commGalGrid_t *thisCommGalGrid, grid_t *thisGrid);
void resort_gal_zreion_photHI(commGalGrid_t *thisCommGalGrid, int32_t *indexArray);

/* ------------------------------------------------------------------------*/
/* FUNCTIONS TO STORE GRID PROPERTIES IN GALAXIES OF NEXT SNAPSHOT */
/* ------------------------------------------------------------------------*/

void put_grid_properties_to_galaxies_next_snapshot(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, commGalGrid_t *thisCommGalGrid);

#endif