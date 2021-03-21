#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "outgal.h"
#include "phys_const.h"
#include "zas_lists.h"

int get_num_snapshots(outgtree_t **thisTreeList)
{
  /* getting snapnumbers and scalefactors */
  int numSnaps = thisTreeList[0]->galaxies[thisTreeList[0]->numGal-1].snapnumber + 1;
#ifdef MPI
  int recvNumSnaps = 0;
  MPI_Allreduce(&numSnaps, &recvNumSnaps, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  
  numSnaps = recvNumSnaps;
#endif
  
  return numSnaps;
}

float *get_scalefactors(outgtree_t **thisTreeList, int numTrees)
{
  outgtree_t *thisTree = NULL;
  outgal_t *thisGal = NULL;
  int numSnaps = get_num_snapshots(thisTreeList);
  float *scalefactor = allocate_array_float(numSnaps, "scalefactor");
  
  for(int tree=0; tree<numTrees; tree++)
  {
    thisTree = thisTreeList[tree];
    for(int gal=0; gal<thisTree->numGal; gal++)
    {
      thisGal = &(thisTree->galaxies[gal]);
      scalefactor[thisGal->snapnumber] = thisGal->scalefactor;
    }
    if(scalefactor[0] > 0.)
    {
      break;
    }
  }
  
#ifdef MPI
  float *recvScalefactor = allocate_array_float(numSnaps, "recvScalefactor");
  MPI_Allreduce(scalefactor, recvScalefactor, numSnaps, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  
  for(int snap=0; snap<numSnaps; snap++)
    scalefactor[snap] = recvScalefactor[snap];
  
  free(recvScalefactor);
#endif

  return scalefactor;
}

float *get_redshifts(float *scalefactors, int numSnaps)
{
  float *redshifts = allocate_array_float(numSnaps, "times");
  
  for(int i=0; i<numSnaps; i++)
    redshifts[i] = 1./scalefactors[i] - 1.;
  
  return redshifts;
}

float *get_times_from_redshifts(float *redshifts, int numSnaps, double h, double omega_m, double omega_l)
{
  float *times = allocate_array_float(numSnaps, "times");

  for(int i=0; i<numSnaps; i++)
    times[i] = calc_time_from_redshift((double)(redshifts[i]), 1.e10, h, omega_m, omega_l);
  
  return times;
}

float calc_time_from_redshift(double zmin, double zmax, double h, double omega_m, double omega_l)
{
  double prefactor = 2.*Mpc_cm/(3*h*1.e7*sqrt(omega_l));
  double tmp = sqrt(omega_l/omega_m);
  
  return prefactor*(asinh(tmp*pow(1.+zmin, -1.5)) - asinh(tmp*pow(1.+zmax, -1.5)));
}

int find_snapshoft_from_redshift(float redshift, int numSnaps, float *redshifts)
{
  int snap = 0;
  float smallerValue = 0., largerValue = 0.;
  
  /* check border values */
  if(redshift >= redshifts[0])
    snap = 0;
  else if(redshift<= redshifts[numSnaps-1])
    snap = numSnaps-1;
  else
  {
    for(int i=0; i<numSnaps-1; i++)
    {
      largerValue = redshifts[i];
      smallerValue = redshifts[i+1];
      if(smallerValue <= redshift && redshift <= largerValue)
      {
        if(fabs(redshift-smallerValue) < fabs(largerValue-redshift))
          snap = i+1;
        else
          snap = i;
        break;
      }
    }
  }
  
  return snap;
}