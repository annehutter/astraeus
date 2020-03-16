#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "domain.h"

domain_t *initDomain(int nbins, int thisRank, int size)
{
#ifdef MPI
    ptrdiff_t alloc_local, local_n0, local_0_start;
    int tmp_local_n0;
#else
    ptrdiff_t local_n0, local_0_start;
#endif
    
    domain_t *newDomain;
    newDomain = malloc(sizeof(domain_t));
    if(newDomain == NULL)
    {
        fprintf(stderr, "ERROR: initDomain: Not enough memory to allocate domain.\n");
        exit(EXIT_FAILURE);
    }
    
#ifdef MPI
    fftw_mpi_init();
    alloc_local = fftw_mpi_local_size_3d(nbins, nbins, nbins, MPI_COMM_WORLD, &local_n0, &local_0_start);
#else
    local_n0 = nbins;
    local_0_start = 0;
#endif
    
    newDomain->nbins = nbins;
    newDomain->local_n0 = local_n0;
    newDomain->local_0_start = local_0_start;
    
#ifdef MPI
    if(size == 2)
    {
        newDomain->numNeighbours = 1;
        newDomain->listNeighbourRanks = (int*) malloc(newDomain->numNeighbours * sizeof(int));
        if(newDomain->listNeighbourRanks == NULL)
        {
            fprintf(stderr, "ERROR: initDomain: Not enough memory to allocate listNeighbourRanks.\n");
            exit(EXIT_FAILURE);
        }
        if(thisRank == 0) newDomain->listNeighbourRanks[0] = 1;
        else if(thisRank == 1) newDomain->listNeighbourRanks[0] = 0;
        else exit(EXIT_FAILURE);
    }
    else
    {
        newDomain->numNeighbours = 2;
        newDomain->listNeighbourRanks = (int*) malloc(newDomain->numNeighbours * sizeof(int));
        if(newDomain->listNeighbourRanks == NULL)
        {
            fprintf(stderr, "ERROR: initDomain: Not enough memory to allocate listNeighbourRanks.\n");
            exit(EXIT_FAILURE);
        }
        newDomain->listNeighbourRanks[0] = (thisRank - 1 + size)%size;
        newDomain->listNeighbourRanks[1] = (thisRank + 1 + size)%size;
    }
#else
    newDomain->numNeighbours = 0;
    newDomain->listNeighbourRanks = NULL;
#endif
    
    newDomain->thisRank = thisRank;
    newDomain->size = size;
    
    
    newDomain->local_n0_list = (int*) malloc(newDomain->size * sizeof(int));
#ifdef MPI
    for(int i=0; i<newDomain->size; i++)
    {
        if(thisRank == i)
        {
            tmp_local_n0 = local_n0;
        }
        MPI_Bcast(&tmp_local_n0, 1, MPI_INT, i, MPI_COMM_WORLD);
        newDomain->local_n0_list[i] = tmp_local_n0;
    }
#else
    newDomain->local_n0_list[0] = thisRank;
#endif

#ifdef MPI
    fftw_mpi_cleanup();
#endif
    
    return newDomain;
}


void deallocate_domain(domain_t *thisDomain)
{
    if(thisDomain->listNeighbourRanks != NULL) free(thisDomain->listNeighbourRanks);
    if(thisDomain->local_n0_list != NULL) free(thisDomain->local_n0_list);
    
    free(thisDomain);
}

int is_cell_in_domain(domain_t *thisDomain, int z_int)  // z_int is number of cells in z-direction
{
    int lower_limit = thisDomain->local_0_start;
    int upper_limit = thisDomain->local_0_start + thisDomain->local_n0;
    
    if(z_int >= lower_limit && z_int < upper_limit)
    {
        return 1;
    }else{
        return 0;
    }
}

int is_pos_in_domain(domain_t *thisDomain, double z)    // z is relative location in box along z-direction (z_int/nbins)
{
    double lower_limit = (double)thisDomain->local_0_start / (double)thisDomain->nbins;
    double upper_limit = (double)(thisDomain->local_0_start + thisDomain->local_n0) / (double)thisDomain->nbins;
    
    if(z >= lower_limit && z < upper_limit)
    {
        return 1;
    }else{
        return 0;
    }
}

int calc_rank_from_pos(domain_t *thisDomain, double z)
{
    int z_int = z*thisDomain->nbins;
    int rank = 0;
    int tmp = 0;
    
    for(int i=0; i<thisDomain->size; i++)
    {
        tmp += thisDomain->local_n0_list[i];
        if(z_int < tmp)
        {
            rank = i;
            break;
        }
    }
    
    return rank;
}

#ifdef MPI
int calc_next_domain(domain_t *thisDomain, double z)
#else
int calc_next_domain(domain_t *thisDomain)
#endif
{
#ifdef MPI
    int z_int = (int)(z * thisDomain->nbins);
    
    if(z_int < thisDomain->local_0_start && z_int != 0)
    {
        return (thisDomain->thisRank - 1 + thisDomain->size)%thisDomain->size;
    }else{
        return (thisDomain->thisRank + 1 + thisDomain->size)%thisDomain->size;
    }
#else
    return thisDomain->thisRank;
#endif
}

int is_this_rank_neighbor(domain_t *thisDomain, int head)
{
    int tmp = 0;
    int numNeighbours = thisDomain->numNeighbours;
    int *listNeighbourRanks = thisDomain->listNeighbourRanks;
    
    for(int i=0; i<numNeighbours; i++)
    {
        if(head == listNeighbourRanks[i])
        {
            tmp = 1;
            break;
        }
    }
    return tmp;
}
