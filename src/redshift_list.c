#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include "utils.h"
#include "cifog/phys_const.h"

#include "redshift_list.h"

int get_num_lines(char *fileName)
{
  char line[MAXLENGTH];
  FILE *f = NULL;
  
  int num = 0;
  int snap = 0;
  float redshift = 0.;
  float scalefactor = 0.;
  
  if( (f=fopen(fileName, "r")) == NULL)
  {
    fprintf(stderr, "Could not open file\n");
    exit(EXIT_FAILURE);
  }
  
  while((fgets(line, sizeof(line), f)))
  {
    sscanf(line, "%d\t%f\t%f", &snap, &redshift, &scalefactor);
    num++;
  }
  
  fclose(f);
  
  return num;
}

redshiftlist_t *read_redshift_snap_list(char *fileName)
{
  char line[MAXLENGTH];
  FILE *f = NULL;
  
  int snap = 0;
  double redshift = 0.;
  float scalefactor = 0.;
  double *redshifts = NULL;
  redshiftlist_t *thisRedshiftList = malloc(sizeof(redshiftlist_t));
  if(thisRedshiftList == NULL)
  {
    fprintf(stderr, "Could not allocate redshiftlist\n");
    exit(EXIT_FAILURE);
  }
    
  int numLines = get_num_lines(fileName);
  
  if( (f=fopen(fileName, "r")) == NULL)
  {
    fprintf(stderr, "Could not open file\n");
    exit(EXIT_FAILURE);
  }
  
  /* SKIP HEADER */
  for(int i=0; i<2; i++)
    fgets(line, sizeof(line), f);
  
  redshifts = allocate_array_double(numLines - 2, "redshiftList");
  int counter = 0;
  
  while((fgets(line, sizeof(line), f)))
  {
    sscanf(line, "%d\t%lf\t%f", &snap, &redshift, &scalefactor);
    redshifts[counter] = redshift;
    counter++;
  }
     
  fclose(f);
  
  thisRedshiftList->numRedshifts = counter;
  thisRedshiftList->redshifts = redshifts;
  
  return thisRedshiftList;
}

void deallocate_redshiftlist(redshiftlist_t *thisRedshiftList)
{
  if(thisRedshiftList->redshifts != NULL) free(thisRedshiftList->redshifts);
  free(thisRedshiftList);
}

double calc_deltaRedshift(redshiftlist_t *thisRedshiftList, int snap)
{
  double result = 0.;
  
  if(snap + 1 < thisRedshiftList->numRedshifts)
    result = thisRedshiftList->redshifts[snap] - thisRedshiftList->redshifts[snap+1];
  
  return result;
}

int get_snap_from_redshift(redshiftlist_t *thisRedshiftList, double redshift)
{
  int snap = thisRedshiftList->redshifts[thisRedshiftList->numRedshifts - 1];
  
  for(int i=0; i<thisRedshiftList->numRedshifts; i++)
  {
    if(thisRedshiftList->redshifts[i] < redshift)
    {
      if(i > 0)
      {
        double diffLow = redshift - thisRedshiftList->redshifts[i];
        double diffUp = thisRedshiftList->redshifts[i-1] - redshift;
        if(diffUp < diffLow)
          snap = i-1;
        else
          snap = i;
        break;
      }
      else
      {
        snap = 0;
        break;
      }
    }
  }
  
  if(snap < 0) 
    snap = 0;
  
  return snap;
}

float *get_times_from_redshifts(redshiftlist_t *thisRedshiftList, int num, double h, double omega_m, double omega_l)
{
  assert(num <= thisRedshiftList->numRedshifts);
  float *times = allocate_array_float(num, "times");

  for(int i=0; i<num; i++)
  {
    times[i] = calc_time_from_redshift(thisRedshiftList->redshifts[i], 1.e10, h, omega_m, omega_l);
//     printf("z = %e\t t = %e\n", thisRedshiftList->redshifts[i], times[i]/(3600.*24.*365.*1.e6));
  }
  
  return times;
}

float calc_time_from_redshift(double zmin, double zmax, double h, double omega_m, double omega_l)
{
        double prefactor = 2.*Mpc_cm/(3*h*1.e7*sqrt(omega_l));
        double tmp = sqrt(omega_l/omega_m);
        
        return prefactor*(asinh(tmp*pow(1.+zmin, -1.5)) - asinh(tmp*pow(1.+zmax, -1.5)));
}