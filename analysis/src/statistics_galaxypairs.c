#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "statistics_galaxypairs.h"

void send_recv_indeces_pos_primaryGalaxies(int *numPrimaryGalOnRanks, int **primaryGalIndex, double **primaryGalPosx, double **primaryGalPosy, double **primaryGalPosz, int thisRank, int size)
{
  MPI_Status status;
  
  /* first send how many primary galaxies there are in total */
  int totNumPrimaryGalOnRanks = 0;
  int *cumNumPrimaryGalOnRanks = allocate_array_int(size, "cumNumPrimaryGalOnRanks");
  for(int rank=0; rank<size; rank++)
  {
    totNumPrimaryGalOnRanks += numPrimaryGalOnRanks[rank];
    if(rank == 0)
      cumNumPrimaryGalOnRanks[rank] = 0;
    else
      cumNumPrimaryGalOnRanks[rank] = cumNumPrimaryGalOnRanks[rank-1] + numPrimaryGalOnRanks[rank-1];
  }
  
  /* create arrays for primary galaxies on each rank */
  int *totPrimaryGalIndex = allocate_array_int(totNumPrimaryGalOnRanks, "totPrimaryGalIndex");
  double *totPrimaryGalPosx = allocate_array_double(totNumPrimaryGalOnRanks, "totPrimaryGalPosx");
  double *totPrimaryGalPosy = allocate_array_double(totNumPrimaryGalOnRanks, "totPrimaryGalPosy");
  double *totPrimaryGalPosz = allocate_array_double(totNumPrimaryGalOnRanks, "totPrimaryGalPosz");
  
  /* broadcast primary galaxies from all ranks to all ranks */
  for(int head=0; head<size; head++)
  {
    for(int destRank=0; destRank<size; destRank++)
    {
      if(destRank != head && thisRank == head)
      {
        if(numPrimaryGalOnRanks[thisRank] > 0)
        {
          MPI_Send(*primaryGalIndex, numPrimaryGalOnRanks[thisRank], MPI_INT, destRank, 100, MPI_COMM_WORLD);
          MPI_Send(*primaryGalPosx, numPrimaryGalOnRanks[thisRank], MPI_DOUBLE, destRank, 101, MPI_COMM_WORLD);
          MPI_Send(*primaryGalPosy, numPrimaryGalOnRanks[thisRank], MPI_DOUBLE, destRank, 102, MPI_COMM_WORLD);
          MPI_Send(*primaryGalPosz, numPrimaryGalOnRanks[thisRank], MPI_DOUBLE, destRank, 103, MPI_COMM_WORLD);
        }
      }
      
      if(destRank == head && destRank == thisRank)
      {
        for(int i=0; i<numPrimaryGalOnRanks[head]; i++)
        {
          totPrimaryGalIndex[cumNumPrimaryGalOnRanks[head] + i] = (*primaryGalIndex)[i];
        }
      }
    }
    if(thisRank != head)
    {
      if(numPrimaryGalOnRanks[head] > 0)
      {
        MPI_Recv(&(totPrimaryGalIndex[cumNumPrimaryGalOnRanks[head]]), numPrimaryGalOnRanks[head], MPI_INT, head, 100, MPI_COMM_WORLD, &status);
        MPI_Recv(&(totPrimaryGalPosx[cumNumPrimaryGalOnRanks[head]]), numPrimaryGalOnRanks[head], MPI_DOUBLE, head, 101, MPI_COMM_WORLD, &status);
        MPI_Recv(&(totPrimaryGalPosy[cumNumPrimaryGalOnRanks[head]]), numPrimaryGalOnRanks[head], MPI_DOUBLE, head, 102, MPI_COMM_WORLD, &status);
        MPI_Recv(&(totPrimaryGalPosz[cumNumPrimaryGalOnRanks[head]]), numPrimaryGalOnRanks[head], MPI_DOUBLE, head, 103, MPI_COMM_WORLD, &status);
      }
    }
  }
  
  free(*primaryGalIndex);
  free(*primaryGalPosx);
  free(*primaryGalPosy);
  free(*primaryGalPosz);
  free(cumNumPrimaryGalOnRanks);
  
  *primaryGalIndex = totPrimaryGalIndex;
  *primaryGalPosx = totPrimaryGalPosx;
  *primaryGalPosy = totPrimaryGalPosy;
  *primaryGalPosz = totPrimaryGalPosz;
}

double get_distanceSqared_galaxy_pair(double posx1, double posy1, double posz1, double posx2, double posy2, double posz2)
{
  double result = (posx1 - posx2)*(posx1 - posx2) + (posy1 - posy2)*(posy1 - posy2) + (posz1 - posz2)*(posz1 - posz2);
  
  return result;
}

void search_secondaryGalaxies(int numPrimaryGal, int *primaryGalOriginRanks, double *primaryGalPosx, double *primaryGalPosy, double *primaryGalPosz, int numGal, double *selectionProperty, double minSelectionProperty, double maxSelectionProperty, double maxDistance, double *galPosx, double *galPosy, double *galPosz, int *galIndex, int ***galIndexAroundPrimaryGal, double ***galDistanceFromPrimaryGal, int ***galOriginRankOfPrimaryGal)
{
  double distance = 0.;
  int **galIndexPrimaryGal = allocate_array_int_pointer(numPrimaryGal, "galIndexPrimaryGal");
  double **galDistancePrimaryGal = allocate_array_double_pointer(numPrimaryGal, "galDistancePrimaryGal");
  int **galOriginRankPrimaryGal = allocate_array_int_pointer(numPrimaryGal, "galOriginRankPrimaryGal");
  int *numAllocated = allocate_array_int(numPrimaryGal, "numAllocated");
  for(int primaryGal=0; primaryGal<numPrimaryGal; primaryGal++)
  {
    galIndexPrimaryGal[primaryGal] = allocate_array_int(MAXLENGTH, "galIndexPrimaryGal[primaryGal]");
    galIndexPrimaryGal[primaryGal][0] = 0;
    galDistancePrimaryGal[primaryGal] = allocate_array_double(MAXLENGTH, "galDistancePrimaryGal[primaryGal]");
    galDistancePrimaryGal[primaryGal][0] = 0.;
    galOriginRankPrimaryGal[primaryGal] = allocate_array_int(MAXLENGTH, "galOriginRankPrimaryGal[primaryGal]");
    galOriginRankPrimaryGal[primaryGal][0] = 0;
    numAllocated[primaryGal] = MAXLENGTH;
  }
  
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  for(int gal=0; gal<numGal; gal++)
  {
    if(selectionProperty[gal] >= minSelectionProperty && selectionProperty[gal] < maxSelectionProperty)
    {
      for(int primaryGal=0; primaryGal<numPrimaryGal; primaryGal++)
      {
        distance = get_distanceSqared_galaxy_pair(primaryGalPosx[primaryGal], primaryGalPosy[primaryGal], primaryGalPosz[primaryGal], galPosx[gal], galPosy[gal], galPosz[gal]);
        if(distance <= maxDistance*maxDistance)
        {
          galIndexPrimaryGal[primaryGal][0]++;
          galIndexPrimaryGal[primaryGal][galIndexPrimaryGal[primaryGal][0]] = galIndex[gal];
          galDistancePrimaryGal[primaryGal][0] += 1.;
          galDistancePrimaryGal[primaryGal][galIndexPrimaryGal[primaryGal][0]] = sqrt(distance);
          galOriginRankPrimaryGal[primaryGal][0]++;
          galOriginRankPrimaryGal[primaryGal][galIndexPrimaryGal[primaryGal][0]] = primaryGalOriginRanks[primaryGal];
          
          if(galIndexPrimaryGal[primaryGal][0] >= numAllocated[primaryGal])
          {
            numAllocated[primaryGal] = numAllocated[primaryGal] * 1.5;
            galIndexPrimaryGal[primaryGal] = realloc(galIndexPrimaryGal[primaryGal], numAllocated[primaryGal] * sizeof(int));
            galDistancePrimaryGal[primaryGal] = realloc(galDistancePrimaryGal[primaryGal], numAllocated[primaryGal] * sizeof(double));
            galOriginRankPrimaryGal[primaryGal] = realloc(galOriginRankPrimaryGal[primaryGal], numAllocated[primaryGal] * sizeof(int));
          }
        }
      }
    }
  }
    
  for(int primaryGal=0; primaryGal<numPrimaryGal; primaryGal++)
  {
    galIndexPrimaryGal[primaryGal] = realloc(galIndexPrimaryGal[primaryGal], (galIndexPrimaryGal[primaryGal][0] + 1) * sizeof(int));
    galDistancePrimaryGal[primaryGal] = realloc(galDistancePrimaryGal[primaryGal], (galIndexPrimaryGal[primaryGal][0] + 1) * sizeof(double));
    galOriginRankPrimaryGal[primaryGal] = realloc(galOriginRankPrimaryGal[primaryGal], (galOriginRankPrimaryGal[primaryGal][0] + 1) * sizeof(int));
  }
  
#ifdef MPI 
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  free(numAllocated);
  *galIndexAroundPrimaryGal = galIndexPrimaryGal;
  *galDistanceFromPrimaryGal = galDistancePrimaryGal;
  *galOriginRankOfPrimaryGal = galOriginRankPrimaryGal;
}

void write_galaxypairs_to_file(int numPrimaryGal, int *primaryGalIndex, int thisRank, int size, int *numPrimaryGalOnRanks, int numGal, double *Mvir, double *MgasIni, double *Mgas, double *Mstar, double *density, double *XHII, int **galIndexAroundPrimaryGal, double **galDistanceFromPrimaryGal, int **galOriginRankOfPrimaryGal, char *basename)
{
  int galIndex = 0;
  char fileEnding[10];
  sprintf(fileEnding, "_properties_%d.dat", thisRank);
  char *filename = concat_strings(2, basename, fileEnding);
    
  FILE *f = fopen(filename, "w");
  if(f == NULL)
  {
    fprintf(stderr, "Could not open output file %s\n", filename);
    exit(EXIT_FAILURE);
  }
  
  fprintf(f, "# Primary (0) or secondary (1) galaxy\t rank\t index of primary galaxy\tdistance [h^-1 Mpc]\t Mvir [Msun]\t MgasIni [Msun]\t Mgas [Msun]\t Mstar [Msun]\t log10(1+delta)\t XHII\n");
  
  int counter=0;
  for(int rank=0; rank<thisRank; rank++)
  {
    counter += numPrimaryGalOnRanks[rank];
  }
  
  for(int primaryGal=counter; primaryGal<counter+numPrimaryGalOnRanks[thisRank]; primaryGal++)
  {
    galIndex = primaryGalIndex[primaryGal];
    if(galIndex >= numGal) printf("rank %d: galIndex = %d\t numGal = %d\n", thisRank, galIndex, numGal);
    fprintf(f, "%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 0, thisRank, primaryGalIndex[primaryGal], 0.,  Mvir[galIndex], MgasIni[galIndex], Mgas[galIndex], Mstar[galIndex], density[galIndex], XHII[galIndex]);
  }
  
  for(int primaryGal=0; primaryGal<numPrimaryGal; primaryGal++)
  {
    for(int gal=0; gal<galIndexAroundPrimaryGal[primaryGal][0]; gal++)
    {
      galIndex = galIndexAroundPrimaryGal[primaryGal][gal + 1];
      if(galIndex >= numGal) printf("rank %d: galIndex = %d\t numGal = %d\n", thisRank, galIndex, numGal);
      fprintf(f, "%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", 1, galOriginRankOfPrimaryGal[primaryGal][gal + 1], primaryGalIndex[primaryGal], galDistanceFromPrimaryGal[primaryGal][gal + 1], Mvir[galIndex], MgasIni[galIndex], Mgas[galIndex], Mstar[galIndex], density[galIndex], XHII[galIndex]);
    }
  }
  
  fclose(f);
  free(filename);
}

void write_galaxypairs_SFHs_to_file(int numPrimaryGal, int *primaryGalIndex, int thisRank, int size, int *numPrimaryGalOnRanks, int numGal, int currSnap, float *redshifts, double *propertyWithHistory, int **galIndexAroundPrimaryGal, double **galDistanceFromPrimaryGal, int **galOriginRankOfPrimaryGal, char *basename)
{
  int galIndex = 0;
  char fileEnding[10];
  sprintf(fileEnding, "_histories_%d.dat", thisRank);
  char *filename = concat_strings(2, basename, fileEnding);
  
  FILE *f = fopen(filename, "w");
  if(f == NULL)
  {
    fprintf(stderr, "Could not open output file %s\n", filename);
    exit(EXIT_FAILURE);
  }
  
  fprintf(f, "# Primary (0) or secondary (1) galaxy\t rank\t index of primary galaxy\t distance []\t assembly history: ");
  for(int snap=0; snap<currSnap+1; snap++)
  {
    fprintf(f, "z = %f\t", redshifts[snap]);
  }
  fprintf(f, "\n");
  
  fprintf(f, "%d\t%d\t%d\t%e\t", 0, 0, 0, 0.);
  for(int snap=0; snap<currSnap+1; snap++)
  {
    fprintf(f, "%f\t", redshifts[snap]);
  }
  fprintf(f, "\n");
  
  int counter=0;
  for(int rank=0; rank<thisRank; rank++)
  {
    counter += numPrimaryGalOnRanks[rank];
  }
  
  for(int primaryGal=counter; primaryGal<counter+numPrimaryGalOnRanks[thisRank]; primaryGal++)
  {
    galIndex = primaryGalIndex[primaryGal];
    if(galIndex >= numGal) printf("rank %d: galIndex = %d\t numGal = %d\n", thisRank, galIndex, numGal);
    fprintf(f, "%d\t%d\t%d\t%e\t", 0, thisRank, primaryGalIndex[primaryGal], 0.);
    for(int snap=0; snap<currSnap+1; snap++)
    {
      fprintf(f, "%e\t", propertyWithHistory[galIndex * (currSnap + 1) + snap]);
    }
    fprintf(f, "\n");
  }
  
  for(int primaryGal=0; primaryGal<numPrimaryGal; primaryGal++)
  {
    for(int gal=0; gal<galIndexAroundPrimaryGal[primaryGal][0]; gal++)
    {
      galIndex = galIndexAroundPrimaryGal[primaryGal][gal + 1];
      if(galIndex >= numGal) printf("rank %d: galIndex = %d\t numGal = %d\n", thisRank, galIndex, numGal);
      fprintf(f, "%d\t%d\t%d\t%e\t", 1, galOriginRankOfPrimaryGal[primaryGal][gal + 1], primaryGalIndex[primaryGal], galDistanceFromPrimaryGal[primaryGal][gal]);
      for(int snap=0; snap<currSnap+1; snap++)
      {
        fprintf(f, "%e\t", propertyWithHistory[galIndex * (currSnap + 1) + snap]);
      }
      fprintf(f, "\n");
    }
  }
  
  fclose(f);
  free(filename);
}

void statistics_on_galaxypairs(int numGal, double *selectionProperty, double minSelectionProperty, double maxSelectionProperty, double *selectionProperty2, double minSelectionProperty2, double maxSelectionProperty2, double maxDistance, int numSnaps, float *redshifts, int currSnap, double *posxProperty, double *posyProperty, double *poszProperty, double *Mvir, double *MgasIni, double *Mgas, double *Mstar, double *density, double *XHII, double *propertyWithHistory, int thisRank, int size, char *basename)
{
  int *galIndex = allocate_array_int(numGal, "galIndex");
  
  /* create list of first selected galaxy sample */
  int *primaryGalIndex = allocate_array_int(MAXLENGTH, "galIndexList");
  double *primaryGalPosx = allocate_array_double(MAXLENGTH, "primaryGalPosx");
  double *primaryGalPosy = allocate_array_double(MAXLENGTH, "primaryGalPosy");
  double *primaryGalPosz = allocate_array_double(MAXLENGTH, "primaryGalPosz");

  int numPrimaryGalOnThisRank = 0; /* number of primary galaxies on this rank */
  int allocatedGal = MAXLENGTH;
  for(int gal=0; gal<numGal; gal++)
  {
    galIndex[gal] = gal;
    if(selectionProperty[gal] >= minSelectionProperty && selectionProperty[gal] < maxSelectionProperty)
    {
      primaryGalIndex[numPrimaryGalOnThisRank] = gal;
      primaryGalPosx[numPrimaryGalOnThisRank] = posxProperty[gal];
      primaryGalPosy[numPrimaryGalOnThisRank] = posyProperty[gal];
      primaryGalPosz[numPrimaryGalOnThisRank] = poszProperty[gal];

      numPrimaryGalOnThisRank++;
      
      if(numPrimaryGalOnThisRank >= allocatedGal)
      {
        primaryGalIndex = realloc(primaryGalIndex, allocatedGal * 1.5 * sizeof(int));
        primaryGalPosx = realloc(primaryGalPosx, allocatedGal * 1.5 * sizeof(double));
        primaryGalPosy = realloc(primaryGalPosy, allocatedGal * 1.5 * sizeof(double));
        primaryGalPosz = realloc(primaryGalPosz, allocatedGal * 1.5 * sizeof(double));
      }
    }
  }
  primaryGalIndex = realloc(primaryGalIndex, numPrimaryGalOnThisRank * sizeof(int));
  primaryGalPosx = realloc(primaryGalPosx, numPrimaryGalOnThisRank * sizeof(double));
  primaryGalPosy = realloc(primaryGalPosy, numPrimaryGalOnThisRank * sizeof(double));
  primaryGalPosz = realloc(primaryGalPosz, numPrimaryGalOnThisRank * sizeof(double));

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
    
  int *numPrimaryGalOnRanks = allocate_array_int(size, "numPimaryGalOnRanks");
#ifdef MPI
  /* communicate first selected galaxy sample to all processors */
  MPI_Allgather(&numPrimaryGalOnThisRank, 1, MPI_INT, numPrimaryGalOnRanks, 1, MPI_INT, MPI_COMM_WORLD);
  
  /* Send galaxies indeces and positions to all ranks */
  send_recv_indeces_pos_primaryGalaxies(numPrimaryGalOnRanks, &primaryGalIndex, &primaryGalPosx, &primaryGalPosy, &primaryGalPosz, thisRank, size);
  
  MPI_Barrier(MPI_COMM_WORLD);
#else
  numPrimaryGalOnRanks[0] = numPrimaryGalOnThisRank;
#endif
  
  /* get secondary galaxies on this rank */
  int numPrimaryGal = 0; // needs to be adjusted!
  for(int rank=0; rank<size; rank++)
  {
    numPrimaryGal += numPrimaryGalOnRanks[rank];
  }
  int *primaryGalOriginRank = allocate_array_int(numPrimaryGal, "primaryGalOriginRank");
  int counter = 0;
  for(int rank=0; rank<size; rank++)
  {
    for(int i=0; i<numPrimaryGalOnRanks[rank]; i++)
    {
      primaryGalOriginRank[counter + i] = rank;
    }
    counter += numPrimaryGalOnRanks[rank];
  } 
    
  /* search second galaxies */
  int **galIndexAroundPrimaryGal = NULL;
  double **galDistanceFromPrimaryGal = NULL;
  int **galOriginRankOfPrimaryGal = NULL;
  
  search_secondaryGalaxies(numPrimaryGal, primaryGalOriginRank, primaryGalPosx, primaryGalPosy, primaryGalPosz, numGal, selectionProperty2, minSelectionProperty2, maxSelectionProperty2, maxDistance, posxProperty, posyProperty, poszProperty, galIndex, &galIndexAroundPrimaryGal, &galDistanceFromPrimaryGal, &galOriginRankOfPrimaryGal);

  /* write out galaxies part of pairs that are located on this rank */
  write_galaxypairs_to_file(numPrimaryGal, primaryGalIndex, thisRank, size, numPrimaryGalOnRanks, numGal, Mvir, MgasIni, Mgas, Mstar, density, XHII, galIndexAroundPrimaryGal, galDistanceFromPrimaryGal, galOriginRankOfPrimaryGal, basename);
  
  write_galaxypairs_SFHs_to_file(numPrimaryGal, primaryGalIndex, thisRank, size, numPrimaryGalOnRanks, numGal, currSnap, redshifts, propertyWithHistory, galIndexAroundPrimaryGal, galDistanceFromPrimaryGal, galOriginRankOfPrimaryGal, basename);
  
  for(int primaryGal=0; primaryGal<numPrimaryGal; primaryGal++)
  {
    free(galIndexAroundPrimaryGal[primaryGal]);
    free(galDistanceFromPrimaryGal[primaryGal]);
    free(galOriginRankOfPrimaryGal[primaryGal]);
  }
  free(galIndexAroundPrimaryGal);
  free(galDistanceFromPrimaryGal);
  free(galOriginRankOfPrimaryGal);
  free(numPrimaryGalOnRanks);
  free(primaryGalOriginRank);
  free(primaryGalIndex);
  free(primaryGalPosx);
  free(primaryGalPosy);
  free(primaryGalPosz);
  free(galIndex);
}