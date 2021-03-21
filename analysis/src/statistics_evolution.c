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
#include "statistics.h"
#include "statistics_evolution.h"

int get_numBins_for_bins(int numGal, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1)
{
  int currSnap1 = 0;
  if(containsHistories1 == 1)
    currSnap1 = currSnap;
    
  /* get minimum and maximum */
  float minValue1 = 0., maxValue1 = 0.;
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty1, &minValue1, &maxValue1, currSnap1);
    
  /* find binning */
  int numBins1 = 0;
  float *values1 = NULL;
  int minIndex1 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  
  free(values1);
  
  return numBins1;
}

int get_minIndex_for_bins(int numGal, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1)
{
  int currSnap1 = 0;
  if(containsHistories1 == 1)
    currSnap1 = currSnap;
    
  /* get minimum and maximum */
  float minValue1 = 0., maxValue1 = 0.;
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty1, &minValue1, &maxValue1, currSnap1);
    
  /* find binning */
  int numBins1 = 0;
  float *values1 = NULL;
  int minIndex1 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  
  free(values1);
  
  return minIndex1;
}

float *get_values_for_bins(int numGal, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1)
{
  int currSnap1 = 0;
  if(containsHistories1 == 1)
    currSnap1 = currSnap;
    
  /* get minimum and maximum */
  float minValue1 = 0., maxValue1 = 0.;
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty1, &minValue1, &maxValue1, currSnap1);
    
  /* find binning */
  int numBins1 = 0;
  float *values1 = NULL;
  int minIndex1 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  
  return values1;
}

int get_numBins_for_fixed_bins(float minValue1, float maxValue1, int binsInLog1, int binsPerMag1)
{   
  /* find binning */
  int numBins1 = 0;
  float *values1 = NULL;
  int minIndex1 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  
  free(values1);
  
  return numBins1;
}

int get_minIndex_for_fixed_bins(float minValue1, float maxValue1, int binsInLog1, int binsPerMag1)
{   
  /* find binning */
  int numBins1 = 0;
  float *values1 = NULL;
  int minIndex1 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  
  free(values1);
  
  return minIndex1;
}

float *get_values_for_fixed_bins(float minValue1, float maxValue1, int binsInLog1, int binsPerMag1)
{   
  /* find binning */
  int numBins1 = 0;
  float *values1 = NULL;
  int minIndex1 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  
  return values1;
}

void calc_1D_histogram_evolution(int numGal, double *property, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int minIndex1, int containsHistories1, int numBins1, int numSnaps, float **histogram, int **histogramNum)
{
  int currSnap1 = 0;
  if(containsHistories1 == 1)
    currSnap1 = currSnap;
  
  float value1 = 0.;
  int valueInt1 = 0;
    
  for(int gal=0; gal<numGal; gal++)
  {
    if(binsInLog1 == 1)
    {
      if(binProperty1[gal * (currSnap1 + 1) + currSnap1] <= 0.)
        value1 = (float)minIndex1 / (float)binsPerMag1;
      else
        value1 = log10(binProperty1[gal * (currSnap1 + 1) + currSnap1]);
    }
    else
      value1 = binProperty1[gal * (currSnap1 + 1) + currSnap1];
    
    valueInt1 = floor(value1 * (float)binsPerMag1 - minIndex1);

    if(valueInt1 < 0 && binsInLog1 == 1)
      valueInt1 = 0;

    assert(valueInt1 >= 0);
    assert(valueInt1 < numBins1);
    
    if(property[gal * (currSnap1 + 1) + currSnap1] >= 0.)
    {
      (*histogram)[valueInt1 * numSnaps + currSnap] += property[gal * (currSnap1 + 1) + currSnap1];
      (*histogramNum)[valueInt1 * numSnaps + currSnap] += 1;
    }
  }
}

void write_1D_histogram_evolution(int numBins, float *values, int endSnap, float *redshifts, float *histogram, int *numHistogram, int thisRank, char *filename)
{
  FILE *f;
  
  if(thisRank == 0)
  {
    f = fopen(filename, "w");
    if(f == NULL)
    {
      fprintf(stderr, "Could not open output file %s\n", filename);
      exit(EXIT_FAILURE);
    }
    for(int valueInt=0; valueInt<numBins; valueInt++)
    {
      for(int snap=0; snap<endSnap; snap++)
      {
        fprintf(f, "%e\t%e\t%e\t%d\n", values[valueInt], redshifts[snap], histogram[valueInt * endSnap + snap], numHistogram[valueInt * endSnap + snap]);
      }
    }
    
    fclose(f);
  }
}