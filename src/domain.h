#ifndef DOMAIN_H
#define DOMAIN_H

typedef struct
{
    int nbins;
    int local_n0;
    int local_0_start;
    int *local_n0_list;
    
    int numNeighbours;
    int *listNeighbourRanks;
    int thisRank;
    int size;
} domain_t;

domain_t *initDomain(int nbins, int thisRank, int size);
void deallocate_domain(domain_t *thisDomain);

int is_cell_in_domain(domain_t *thisDomain, int z_int);  // z_int is number of cells in z-direction
int is_pos_in_domain(domain_t *thisDomain, double z);    // z is relative location in box along z-direction (z_int/nbins)
int calc_rank_from_pos(domain_t *thisDomain, double z);

#ifdef MPI
int calc_next_domain(domain_t *thisDomain, double z);
#else
int calc_next_domain(domain_t *thisDomain);
#endif

int is_this_rank_neighbor(domain_t *thisDomain, int head);

#endif
