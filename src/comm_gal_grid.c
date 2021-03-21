#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>

#ifdef MPI
#include <fftw3-mpi.h>
#include <mpi.h>
#else
#include <fftw3.h>
#endif

#include "utils.h"
#include "dconfObj.h"
#include "gal_gtree.h"
#include "domain.h"

#include "cifog/confObj.h"
#include "cifog/grid.h"

#include "sort.h"
#include "comm_gal_grid_struct.h"
#include "comm_gal_grid.h"

/* ------------------------------------------------------------------------*/
/* MAIN FUNCTION TO GET ZREION & PHOTHI FOR GALAXIES IN NEXT SNAPSHOT      */
/* ------------------------------------------------------------------------*/

void map_grid_to_gal(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, dconfObj_t simParam, domain_t *thisDomain, grid_t *thisGrid)
{
  /* ----------------------------------------------------------*/
  /* get positions and grid ranks of galaxies in next snapshot */
  /* ----------------------------------------------------------*/
  
  /* initialize communication arrays for gal-grid communication (allocation happens within function) */
  commGalGrid_t *thisCommGalGrid = initCommGalGrid_pos(numGtrees*10);
  int32_t *gridRanks = allocate_array_int32_t(numGtrees*10, "gridRanks");
  
  /* get positions and grid ranks of galaxies */
  get_galaxies_next_snapshot(numGtrees, thisGtreeList, minSnap, maxSnap, simParam, thisDomain, thisCommGalGrid, &gridRanks);
  
  /* ----------------------------------------------------------*/
  /* get grid properties at galaxy positions  */
  /* ----------------------------------------------------------*/ 
  get_grid_properties_galaxies_next_snapshot(thisDomain, thisGrid, thisCommGalGrid, &gridRanks);

  /* ----------------------------------------------------------*/
  /* store grid properties in galaxies (galaxy struct)  */
  /* ----------------------------------------------------------*/ 
  put_grid_properties_to_galaxies_next_snapshot(numGtrees, thisGtreeList, minSnap, maxSnap, thisCommGalGrid);

  /* ----------------------------------------------------------*/
  /* deallocation  */
  /* ----------------------------------------------------------*/ 
  
  deallocate_commGalGrid(thisCommGalGrid);
  if(gridRanks != NULL) free(gridRanks);
}


/* ------------------------------------------------------------------------*/
/* FUNCTIONS TO GET POSITIONS AND GRID RANKS FOR GALAXIES IN NEXT SNAPSHOT */
/* ------------------------------------------------------------------------*/

int32_t get_gal_position_index(gal_t *thisGal, dconfObj_t simParam, domain_t *thisDomain)
{
  int32_t nbins = thisDomain->nbins;
  double inv_boxsize = simParam->inv_boxsize;
  int32_t x = 0, y = 0, z = 0;

  x = thisGal->pos[0] * inv_boxsize * nbins;  // assumes pos[i] to arange between 0 and 1
  y = thisGal->pos[1] * inv_boxsize * nbins;
  z = thisGal->pos[2] * inv_boxsize * nbins;
  if(x >= nbins)
    x = nbins - 1;
  if(y >= nbins)
    y = nbins - 1;
  if(z >= nbins)
    z = nbins - 1;
  
  return z*nbins*nbins + y*nbins + x;
}

int32_t get_gal_gridRank(gal_t *thisGal, dconfObj_t simParam, domain_t *thisDomain)
{
  double inv_boxsize = simParam->inv_boxsize;

  return calc_rank_from_pos(thisDomain, thisGal->pos[2]*inv_boxsize);
}

/* function to get the position of galaxies in next snapshot */
void get_galaxies_next_snapshot(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, dconfObj_t simParam, domain_t *thisDomain, commGalGrid_t *thisCommGalGrid, int32_t **gridRanks)
{
  gtree_t **theseGtrees = *thisGtreeList;
  gtree_t *thisGtree = NULL;
  gal_t *thisGal;
  int32_t gal = 0;
  int32_t snapGal = 0;
  int32_t numGal = 0;
    
  int32_t numGalInSnap = 0;
  int numGalRealloc = 0;
  
  /* LOOP OVER TREES */
  for(int gtree=0; gtree<numGtrees; gtree++)
  {
    thisGtree = theseGtrees[gtree];
    numGal = thisGtree->numGal;
    gal = thisGtree->walker;
    
    if(gal < numGal) 
    {
      snapGal = thisGtree->galaxies[gal].snapnumber;

      /* LOOP OVER SNAPSHOTS */
      for(int snap=minSnap+1; snap<maxSnap+1; snap++)
      {
        if(snapGal == snap)
        {
          /* LOOP OVER GALAXIES */
          thisGal = &(thisGtree->galaxies[gal]);
          while(thisGal->snapnumber == snap)
          { 
            thisCommGalGrid->pos[thisCommGalGrid->numGalWritten] = get_gal_position_index(thisGal, simParam, thisDomain);
            (*gridRanks)[thisCommGalGrid->numGalWritten] = get_gal_gridRank(thisGal, simParam, thisDomain);
            thisCommGalGrid->numGalWritten = thisCommGalGrid->numGalWritten + 1;
            
            if(thisCommGalGrid->numGalWritten >= thisCommGalGrid->numGal)
            {
              if(thisCommGalGrid->numGal > 3)
                numGalRealloc = thisCommGalGrid->numGal * 1.3;
              else
                numGalRealloc = thisCommGalGrid->numGal + 1;
              reallocCommGalGrid_pos(&thisCommGalGrid, numGalRealloc);
              *gridRanks = realloc(*gridRanks, sizeof(int32_t)*numGalRealloc);
            }
            
            gal++;
            numGalInSnap++;
            
            /* stop when reaching the end of the tree */
            if(gal == numGal)
              break;
            
            /* go to the next galaxy */
            thisGal = &(thisGtree->galaxies[gal]);
          }
          if(gal < numGal)
            snapGal = thisGtree->galaxies[gal].snapnumber;
        }
      }
    }
  }
  
//   printf("thisCommGalGrid->numGalWritten = %d\n", thisCommGalGrid->numGalWritten);
  
  /* compress arrays to number of galaxies */
  reallocCommGalGrid_pos(&thisCommGalGrid, thisCommGalGrid->numGalWritten);
  *gridRanks = realloc(*gridRanks, sizeof(int32_t) * thisCommGalGrid->numGalWritten);
}

/* ------------------------------------------------------------------------*/
/* FUNCTIONS TO GET GRID PROPERTIES FOR GALAXIES IN NEXT SNAPSHOT */
/* ------------------------------------------------------------------------*/

/* main function to get grid properties for galaxies (including all communications) */
void get_grid_properties_galaxies_next_snapshot(domain_t *thisDomain, grid_t *thisGrid,  commGalGrid_t *thisCommGalGrid, int32_t **gridRanks)
{
  /* numGal: number of galaxies on this rank */
  /* numGalGridRanks: number of galaxies on grid ranks*/
  /* numGalGalRanks: number of galaxies received from each grid rank (adds up to numGal) */
  
  int numGal = thisCommGalGrid->numGal;
  
  /* create index array */
  int32_t *indexArray = create_index_array(numGal);
  
  /* sort galaxies according to their grid ranks (preparation for sending) */
  quicksort(*gridRanks, 0, numGal-1, indexArray);
  sort_gal_positions(thisCommGalGrid, indexArray);

  /* get how many galaxies are on grid ranks and allocate space for galaxy ranks */
  int32_t *numGal_gridRanks = get_numGal_gridRanks(numGal, *gridRanks, thisDomain->size);
  int32_t *numGal_galRanks = allocate_array_int32_t(thisDomain->size, "numGal_galRanks");
  
  /* communicate galaxy positions to grid ranks */
  int32_t *toRecvPos = NULL;
  send_recv_pos(thisDomain, numGal_gridRanks, &numGal_galRanks, thisCommGalGrid->pos, &toRecvPos);
  if(thisCommGalGrid->pos != NULL) 
  {
    free(thisCommGalGrid->pos);
    thisCommGalGrid->numGal = sum_numGal_ranks(thisDomain->size, numGal_galRanks);
    thisCommGalGrid->numGalWritten = thisCommGalGrid->numGal;
    thisCommGalGrid->pos = toRecvPos;
  }

  /* get grid values at galaxy positions */
  allocCommGalGrid_zreion_photHI(&thisCommGalGrid, thisCommGalGrid->numGal);
  get_gal_zreion_photHI(thisCommGalGrid, thisGrid);
  
  /* send grid values at galaxy positions to galaxy ranks */
  int32_t *recvNumGal_gridRanks = allocate_array_int32_t(thisDomain->size, "recvNumGal_gridRanks");
  float *toRecvZreion = NULL;
  float *toRecvPhotHI = NULL;
  send_recv_zreion_photHI(thisDomain, numGal_galRanks, &recvNumGal_gridRanks, thisCommGalGrid->zreion, thisCommGalGrid->photHI, &toRecvZreion, &toRecvPhotHI);
  check_numGal_gridRanks(thisDomain->size, numGal_gridRanks, recvNumGal_gridRanks);
  if(thisCommGalGrid->zreion != NULL)
  {
    free(thisCommGalGrid->zreion);
    thisCommGalGrid->numGal = sum_numGal_ranks(thisDomain->size, recvNumGal_gridRanks);
    thisCommGalGrid->numGalWritten = thisCommGalGrid->numGal;
    thisCommGalGrid->zreion = toRecvZreion;
  }
  if(thisCommGalGrid->photHI != NULL)
  {
    free(thisCommGalGrid->photHI);
    thisCommGalGrid->numGal = sum_numGal_ranks(thisDomain->size, recvNumGal_gridRanks);
    thisCommGalGrid->numGalWritten = thisCommGalGrid->numGal;
    thisCommGalGrid->photHI = toRecvPhotHI;
  }
  
  /* resort grid values according to index array or figure out smart way of mapping */
  resort_gal_zreion_photHI(thisCommGalGrid, indexArray);

  /* deallocation */
  if(indexArray != NULL) free(indexArray);
  if(numGal_gridRanks != NULL) free(numGal_gridRanks);
  if(numGal_galRanks != NULL) free(numGal_galRanks);
  if(recvNumGal_gridRanks != NULL) free(recvNumGal_gridRanks);
}

int32_t *create_index_array(int numGal)
{
    int32_t *indexArray = allocate_array_int32_t(numGal, "indexArray");

    for(int gal=0; gal<numGal; gal++)
      indexArray[gal] = gal;
    
    return indexArray;
}

void sort_gal_positions(commGalGrid_t *thisCommGalGrid, int32_t *indexArray)
{
  int numGal = thisCommGalGrid->numGalWritten;
  int32_t *posArray = allocate_array_int32_t(numGal, "posArray");
  
  for(int gal=0; gal<numGal; gal++)
    posArray[gal] = thisCommGalGrid->pos[indexArray[gal]];
  
  free(thisCommGalGrid->pos);
  thisCommGalGrid->pos = posArray;
}

/* GET HOW MANY GALAXIES HAVE TO BE SENT TO EACH GRID RANK (works only on a sorted array!) */
int32_t *get_numGal_gridRanks(int numGal, int32_t *gridRanks, int32_t size)
{
  int32_t *numGal_gridRanks = allocate_array_int32_t(size, "numGal_gridRanks");
  
  for(int gal=0; gal<numGal; gal++)
  {
    numGal_gridRanks[gridRanks[gal]]++;
  }
    
  return numGal_gridRanks;
}

/* SEND & RECEIVE POS TO PROCESSORS ACCORDING TO THEIR LOCATION */
#ifdef MPI
void send_recv_pos(domain_t *thisDomain, int32_t *numGal_gridRanks, int32_t **numGal_galRanks, int32_t *toSendPos, int32_t **toRecvPos)
{
  MPI_Status status;
  int size = thisDomain->size;
  int thisRank = thisDomain->thisRank;
  
  int32_t head_offset = 0;
  
  int32_t *recvNumGalOnRank = *numGal_galRanks;
  int32_t counterRecvNumGal = 0;
  
  /* head is the sending rank */
  for(int head=0; head<size; head++)
  {    
    head_offset = 0;
    for(int destRank=0; destRank<size; destRank++)
    {
      if((thisRank == head) && (thisRank != destRank))
      {
        /* send stuff from head to destRank */
        MPI_Send(&numGal_gridRanks[destRank], 1, MPI_INT, destRank, 100, MPI_COMM_WORLD);
        
        if(numGal_gridRanks[destRank] > 0)
        {            
          MPI_Send(&toSendPos[head_offset], numGal_gridRanks[destRank], MPI_INT, destRank, 102, MPI_COMM_WORLD);
        }
      }
      
      if((thisRank == destRank) && (head == destRank))
      {
        /* (send) copy stuff within rank */
        recvNumGalOnRank[head] = numGal_gridRanks[destRank];
        
        counterRecvNumGal += recvNumGalOnRank[head];
        (*toRecvPos) = realloc((*toRecvPos), counterRecvNumGal * sizeof(int32_t));
        
        for(int i=0; i<recvNumGalOnRank[head]; i++)
        {
            (*toRecvPos)[counterRecvNumGal - recvNumGalOnRank[head] + i] = toSendPos[head_offset + i];
        }
      }
      
      head_offset += numGal_gridRanks[destRank];
    }
    
    if(thisRank != head)
    {
      /* receive stuff from head */
      MPI_Recv(&recvNumGalOnRank[head], 1, MPI_INT, head, 100, MPI_COMM_WORLD, &status);
      
      counterRecvNumGal += recvNumGalOnRank[head];
      if(recvNumGalOnRank[head] > 0)
      {
        (*toRecvPos) = realloc((*toRecvPos), counterRecvNumGal * sizeof(int32_t));
              
        MPI_Recv(&(*toRecvPos)[counterRecvNumGal - recvNumGalOnRank[head]], recvNumGalOnRank[head], MPI_INT, head, 102, MPI_COMM_WORLD, &status);
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  numGal_galRanks = &recvNumGalOnRank;
}
#else
void send_recv_pos(domain_t *thisDomain, int32_t *numGal_gridRanks, int32_t **numGal_galRanks, int32_t *toSendPos, int32_t **toRecvPos)
{ 
  *numGal_galRanks = numGal_galRanks;
  *toRecvPos = toSendPos;
}
#endif

int32_t sum_numGal_ranks(int size, int32_t *numGal_ranks)
{
  int32_t sum = 0;
  
  for(int rank=0; rank<size; rank++)
    sum += numGal_ranks[rank];
  
  return sum;
}

/* SEND & RECEIVE POS TO PROCESSORS ACCORDING TO THEIR LOCATION */
#ifdef MPI
void send_recv_zreion_photHI(domain_t *thisDomain, int32_t *numGal_gridRanks, int32_t **numGal_galRanks, float *toSendZreion, float *toSendPhotHI, float **toRecvZreion, float **toRecvPhotHI)
{
  MPI_Status status;
  int size = thisDomain->size;
  int thisRank = thisDomain->thisRank;
  
  int32_t head_offset = 0;
  
  int32_t *recvNumGalOnRank = *numGal_galRanks;
  int32_t counterRecvNumGal = 0;
  
  /* head is the sending rank */
  for(int head=0; head<size; head++)
  {    
    head_offset = 0;
    for(int destRank=0; destRank<size; destRank++)
    {
      if((thisRank == head) && (thisRank != destRank))
      {
        /* send stuff from head to destRank */
        MPI_Send(&numGal_gridRanks[destRank], 1, MPI_INT, destRank, 100, MPI_COMM_WORLD);
        
        if(numGal_gridRanks[destRank] > 0)
        {            
          MPI_Send(&toSendZreion[head_offset], numGal_gridRanks[destRank], MPI_FLOAT, destRank, 101, MPI_COMM_WORLD);
          MPI_Send(&toSendPhotHI[head_offset], numGal_gridRanks[destRank], MPI_FLOAT, destRank, 102, MPI_COMM_WORLD);
        }
      }
      
      if((thisRank == destRank) && (head == destRank))
      {
        /* (send) copy stuff within rank */
        recvNumGalOnRank[head] = numGal_gridRanks[destRank];
        
        counterRecvNumGal += recvNumGalOnRank[head];
        (*toRecvZreion) = realloc((*toRecvZreion), counterRecvNumGal * sizeof(float));
        (*toRecvPhotHI) = realloc((*toRecvPhotHI), counterRecvNumGal * sizeof(float));

        for(int i=0; i<recvNumGalOnRank[head]; i++)
        {
            (*toRecvZreion)[counterRecvNumGal - recvNumGalOnRank[head] + i] = toSendZreion[head_offset + i];
            (*toRecvPhotHI)[counterRecvNumGal - recvNumGalOnRank[head] + i] = toSendPhotHI[head_offset + i];
        }
      }
      
      head_offset += numGal_gridRanks[destRank];
    }
    
    if(thisRank != head)
    {
      /* receive stuff from head */
      MPI_Recv(&recvNumGalOnRank[head], 1, MPI_INT, head, 100, MPI_COMM_WORLD, &status);
      
      counterRecvNumGal += recvNumGalOnRank[head];
      if(recvNumGalOnRank[head] > 0)
      {
        (*toRecvZreion) = realloc((*toRecvZreion), counterRecvNumGal * sizeof(float));
        (*toRecvPhotHI) = realloc((*toRecvPhotHI), counterRecvNumGal * sizeof(float));
        
        MPI_Recv(&(*toRecvZreion)[counterRecvNumGal - recvNumGalOnRank[head]], recvNumGalOnRank[head], MPI_FLOAT, head, 101, MPI_COMM_WORLD, &status);
        MPI_Recv(&(*toRecvPhotHI)[counterRecvNumGal - recvNumGalOnRank[head]], recvNumGalOnRank[head], MPI_FLOAT, head, 102, MPI_COMM_WORLD, &status);
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  numGal_galRanks = &recvNumGalOnRank;
}
#else
void send_recv_zreion_photHI(domain_t *thisDomain, int32_t *numGal_gridRanks, int32_t **numGal_galRanks, float *toSendZreion, float *toSendPhotHI, float **toRecvZreion, float **toRecvPhotHI)
{ 
  *numGal_galRanks = numGal_galRanks;
  *toRecvZreion = toSendZreion;
  *toRecvPhotHI = toSendPhotHI;
}
#endif

int32_t check_numGal_gridRanks(int size, int32_t *numGal_gridRanks, int32_t *recvNumGal_gridRanks)
{
  int result = 1;

  for(int rank=0; rank<size; rank++)
  {  
    if(numGal_gridRanks[rank] != recvNumGal_gridRanks[rank])
      result = 0;
  }

  return result;
}

/* GET GRID VALUES */
void get_gal_zreion_photHI(commGalGrid_t *thisCommGalGrid, grid_t *thisGrid)
{
  int pos = 0;
  int nbins = thisGrid->nbins;
  int local_0_start = thisGrid->local_0_start;
  
  for(int gal=0; gal<thisCommGalGrid->numGal; gal++)
  {
    pos = thisCommGalGrid->pos[gal];
    thisCommGalGrid->zreion[gal] = creal(thisGrid->zreion[pos - local_0_start*nbins*nbins]);
    thisCommGalGrid->photHI[gal] = creal(thisGrid->photHI_zreion[pos - local_0_start*nbins*nbins]);
  }
}

/* RESORT VALUE ARRAYS */
void resort_gal_zreion_photHI(commGalGrid_t *thisCommGalGrid, int32_t *indexArray)
{
  int numGal = thisCommGalGrid->numGal;
  float *zreion = NULL;
  float *photHI = NULL;
  
  if(numGal > 0)
  {
    zreion = allocate_array_float(numGal, "zreion");
    photHI = allocate_array_float(numGal, "photHI");
    
    for(int gal=0; gal<numGal; gal++)
    {
      zreion[indexArray[gal]] = thisCommGalGrid->zreion[gal];
      photHI[indexArray[gal]] = thisCommGalGrid->photHI[gal];
    }
    
    if(thisCommGalGrid->zreion != NULL)
    {
      free(thisCommGalGrid->zreion);
      thisCommGalGrid->zreion = zreion;
    }
    
    if(thisCommGalGrid->photHI != NULL)
    {
      free(thisCommGalGrid->photHI);
      thisCommGalGrid->photHI = photHI;
    }
  }
}

/* ------------------------------------------------------------------------*/
/* FUNCTIONS TO STORE GRID PROPERTIES IN GALAXIES OF NEXT SNAPSHOT */
/* ------------------------------------------------------------------------*/

void put_grid_properties_to_galaxies_next_snapshot(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, commGalGrid_t *thisCommGalGrid)
{
  gtree_t **theseGtrees = *thisGtreeList;
  gtree_t *thisGtree = NULL;
  gal_t *thisGal;
  int32_t gal = 0;
  int32_t snapGal = 0;
  int32_t numGal = 0;
    
  int counter = 0;
  
  /* LOOP OVER TREES */
  for(int gtree=0; gtree<numGtrees; gtree++)
  {
    thisGtree = theseGtrees[gtree];
    numGal = thisGtree->numGal;
    gal = thisGtree->walker;
    
    if(gal < numGal) 
    {
      snapGal = thisGtree->galaxies[gal].snapnumber;

      /* LOOP OVER SNAPSHOTS */
      for(int snap=minSnap+1; snap<maxSnap+1; snap++)
      {
        if(snapGal == snap)
        {
          /* LOOP OVER GALAXIES */
          thisGal = &(thisGtree->galaxies[gal]);
          while(thisGal->snapnumber == snap)
          { 
            if(thisGal->zreion <= 0.) 
            {
              thisGal->zreion = thisCommGalGrid->zreion[counter];
            }
            thisGal->photHI_bg = thisCommGalGrid->photHI[counter];
//             if(gal == thisGtree->walker && gtree == 0) printf("zreion = %e\t photHI = %e\n", thisGal->zreion, thisGal->photHI_bg);
            
            gal++;
            counter++;
            
            /* stop when reaching the end of the tree */
            if(gal == numGal)
              break;
            
            /* go to the next galaxy */
            thisGal = &(thisGtree->galaxies[gal]);
          }
          if(gal < numGal)
            snapGal = thisGtree->galaxies[gal].snapnumber;
        }
      }
    }
  }
}
