#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "dconfObj.h"
#include "cifog/confObj.h"

#include "check_cifogParam.h"

void check_cifogParam(dconfObj_t simParam, confObj_t cifogParam, int thisRank)
{
  if(simParam->startSnap != cifogParam->snapshot_start)
  {
    if(thisRank == 0) printf("Adjusting cifogParam->snapshot_start to %d\n", simParam->startSnap);
    cifogParam->snapshot_start = simParam->startSnap;
  }
  
  if(simParam->gridsize != cifogParam->grid_size)
  {
    if(thisRank == 0) printf("Adjusting cifogParam->grid_size to %d\n", simParam->gridsize);
    cifogParam->grid_size = simParam->gridsize;
  }
  
  if(simParam->boxsize != cifogParam->box_size)
  {
    if(thisRank == 0) printf("Adjusting cifogParam->box_size to %e\n", simParam->boxsize);
    cifogParam->box_size = simParam->boxsize;
  }
  
  if(simParam->hubble_h != cifogParam->h)
  {
    if(thisRank == 0) printf("Adjusting cifogParam->h to %e\n", simParam->hubble_h);
    cifogParam->h = simParam->hubble_h;
  }
  
  if(simParam->omega_b != cifogParam->omega_b)
  {
    if(thisRank == 0) printf("Adjusting cifogParam->omega_b to %e\n", simParam->omega_b);
    cifogParam->omega_b = simParam->omega_b;
  }
  
  if(simParam->omega_m != cifogParam->omega_m)
  {
    if(thisRank == 0) printf("Adjusting cifogParam->omega_m to %e\n", simParam->omega_m);
    cifogParam->omega_m = simParam->omega_m;
  }
  
  if(simParam->omega_l != cifogParam->omega_l)
  {
    if(thisRank == 0) printf("Adjusting cifogParam->omega_l to %e\n", simParam->omega_l);
    cifogParam->omega_l = simParam->omega_l;
  }
}