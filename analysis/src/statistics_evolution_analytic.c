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
#include "statistics_evolution_analytic.h"

void calc_1D_histogram_evolution_analytic(int numGal, double *property, double *numDensCurrSnap, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int minIndex1, int containsHistories1, int numBins1, int numSnaps, float **histogram, int **histogramNum)
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
      (*histogram)[valueInt1 * numSnaps + currSnap] += numDensCurrSnap[gal] * property[gal * (currSnap1 + 1) + currSnap1];
      (*histogramNum)[valueInt1 * numSnaps + currSnap] += numDensCurrSnap[gal];
    }
  }
}
