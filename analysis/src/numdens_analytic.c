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
#include "numdens_analytic.h"

#define MAXHALOS 100000
#define SQR(X) (X)*(X)

char *create_filename_hmf(char *filename, double redshift)
{
  char redshift_name[12];
  if(redshift >= 10.)
    sprintf(redshift_name, "_z%3.1f", redshift);
  else
    sprintf(redshift_name, "_z%2.1f", redshift);
  
  char *newFilename = concat_strings(3, filename, redshift_name, ".txt");
    
  return newFilename;
}

/* function to read respective HMF file */
void read_HMF_file(char *fileName, int *numBins, double **halomass, double **hmf, double hubble_h)
{
  FILE *f = NULL;
  char line[MAXLENGTH];
  int counter = 0;
  double tmp;
  double inv_hubble_h = 1./hubble_h;
  
  *numBins = MAXHALOS;
  *halomass = allocate_array_double(*numBins, "halomass");
  *hmf = allocate_array_double(*numBins, "hmf");
  
  if( (f=fopen(fileName, "r")) == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", fileName);
    exit(EXIT_FAILURE);
  }
    
  while((fgets(line, sizeof(line), f)))
  {
    if(strncmp("#", line, 1) != 0)
    {
      sscanf(line, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &((*halomass)[counter]), &tmp, &tmp, &tmp, &tmp, &tmp, &tmp, &((*hmf)[counter]), &tmp, &tmp, &tmp, &tmp);
      (*halomass)[counter] = log10((*halomass)[counter] * inv_hubble_h);
      (*hmf)[counter] = (*hmf)[counter] * hubble_h * hubble_h * hubble_h;
      counter++;
    }
  }
  
  *numBins = counter;
  *halomass = realloc(*halomass, sizeof(double) * counter);
  *hmf = realloc(*hmf, sizeof(double) * counter);
}


/* interpolate to obtain numDens */
double calc_numDens(int numBins, double *halomass, double *hmf, double thisHalomass, double thisBinwidth)
{
  double result = 0, lowHmf = 0., highHmf = 0.;
  int lowIndex = 0, highIndex = 0;
  double binwidth = 0.;
  
  for(int i=0; i<numBins; i++)
  {
    if(halomass[i] > thisHalomass)
    {
      lowIndex = i-1;
      break;
    }
  }

  for(int i=0; i<numBins; i++)
  {
    if(halomass[i] > (thisHalomass + thisBinwidth))
    {
      highIndex = i-1;
      break;
    }
  }
  
  binwidth = halomass[1] - halomass[0];
  
  if(lowIndex == highIndex)
    result = 0.5 * (hmf[lowIndex+1] + hmf[lowIndex]) * thisBinwidth;
  else
  {     
    lowHmf = hmf[lowIndex] + (hmf[lowIndex+1] - hmf[lowIndex]) / binwidth * (halomass[lowIndex+1] - thisHalomass);
    highHmf = hmf[highIndex] + (hmf[highIndex+1] - hmf[highIndex]) / binwidth * (thisHalomass + thisBinwidth - halomass[highIndex+1]);
    
    result = 0.5 * (lowHmf + hmf[lowIndex+1]) * (halomass[lowIndex+1] - thisHalomass) + 0.5 * (hmf[highIndex] + highHmf) * (thisHalomass + thisBinwidth - halomass[highIndex]);
    
    for(int i=lowIndex+1; i<highIndex; i++)
    {
      result += 0.5 * (hmf[i] + hmf[i+1]) * binwidth;
    }
  }
  
  return result;
}

double *get_numDens(int numGal, double *halomassEndSnap, double thisBinwidth, char *fileName, double hubble_h)
{
  double *numDens = allocate_array_double(numGal, "numDens");
  int numBins = 0;
  double *logHalomass = NULL;
  double *hmf = NULL;
  
  read_HMF_file(fileName, &numBins, &logHalomass, &hmf, hubble_h);
  
  for(int gal=0; gal<numGal; gal++)
  {
    numDens[gal] = calc_numDens(numBins, logHalomass, hmf, log10(halomassEndSnap[gal]), thisBinwidth);
  }
  
  return numDens;
}
  