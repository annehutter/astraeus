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
#include "statistics_analytic.h"

void calc_1D_histogram_analytic(int numGal, double *binProperty1, double *numDensEndSnap, int currSnap, int binsInLog1, int binsPerMag1, int containsHistories1, int cumulative, int thisRank, char *filename)
{
  int currSnap1 = 0;
  if(containsHistories1 == 1)
    currSnap1 = currSnap;
    
  /* get minimum and maximum */
  float minValue1 = 0., maxValue1 = 0.;
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty1, &minValue1, &maxValue1, currSnap1);
    
  /* find binning */
  int numBins1 = 0;
  float value1 = 0;
  float *values1 = NULL;
  int minIndex1 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
    
  /* construct histogram */
  int valueInt1 = 0;
  double *histogramNum = allocate_array_double(numBins1, "histogramNum");
   
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
      
    histogramNum[valueInt1] += numDensEndSnap[gal];
  }
  
  double sumNum = 0.;
  for(int i=0; i<numBins1; i++)
    sumNum += histogramNum[i];
    
  /* communicate histograms */
#ifdef MPI
  double *recvNumHistogram = allocate_array_double(numBins1, "recvNumHistogram");
  
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
      histogramNum[valueInt1] = recvNumHistogram[valueInt1];
  
  free(recvNumHistogram);
  
  double recvSumNum = 0.;
  int recvNumGal = 0;
  
  MPI_Allreduce(&sumNum, &recvSumNum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&numGal, &recvNumGal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  
  numGal = recvNumGal;
  sumNum = recvSumNum;
#endif
  
  double binwidth = 0.;
  if(binsInLog1 == 1)
    binwidth = log10(values1[1]) - log10(values1[0]);
  else
    binwidth = values1[1] - values1[0];
    
  if(cumulative == 1)
  {
    double sum = 0.;
    for(int i=0; i<numBins1; i++)
    {
      sum += histogramNum[i] / binwidth;
      histogramNum[i] = sum;
    }
  }
  else{
    for(int i=0; i<numBins1; i++)
    {
      histogramNum[i] = histogramNum[i] / binwidth;
    }
  }
      
  write_1D_histogram_analytic(numBins1, values1, histogramNum, thisRank, filename);
      
  free(values1);
  free(histogramNum);
}

void calc_2D_histogram_analytic(int numGal, double *binProperty1, double *binProperty2, double *numDensEndSnap, int currSnap, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, char *filename)
{
  int currSnap1 = 0, currSnap2 = 0;
  if(containsHistories1 == 1)
    currSnap1 = currSnap;
  if(containsHistories2 == 1)
    currSnap2 = currSnap;
  
  /* get minimum and maximum */
  float minValue1 = 0., maxValue1 = 0.;
  float minValue2 = 0., maxValue2 = 0.;
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty1, &minValue1, &maxValue1, currSnap1);
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty2, &minValue2, &maxValue2, currSnap2);
  
  /* find binning */
  int numBins1 = 0, numBins2 = 0;
  float value1 = 0, value2 = 0;
  float *values1 = NULL, *values2 = NULL;
  int minIndex1 = 0, minIndex2 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  get_bins(minValue2, maxValue2, binsInLog2, binsPerMag2, &numBins2, &values2, &minIndex2);
  
  if(thisRank == 0)
  {
    printf("2D histogram analytic: minValue1 = %e\t maxValue1 = %e\t numBins1 = %d\n", minValue1, maxValue1, numBins1);
    printf("2D histogram analytic: minValue2 = %e\t maxValue2 = %e\t numBins2 = %d\n", minValue2, maxValue2, numBins2);
  }

  /* construct histogram */
  int valueInt1 = 0, valueInt2 = 0;
  double *histogramNum = allocate_array_double(numBins1 * numBins2, "histogramNum");
  
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
    
    if(binsInLog2 == 1)
    {
      if(binProperty2[gal * (currSnap2 + 1) + currSnap2] <= 0.)
        value2 = (float)minIndex2 / (float)binsPerMag2;
      else
        value2 = log10(binProperty2[gal * (currSnap2 + 1) + currSnap2]);
    }
    else
      value2 = binProperty2[gal * (currSnap2 + 1) + currSnap2];
    
    valueInt1 = floor(value1 * (float)binsPerMag1 - minIndex1);
    valueInt2 = floor(value2 * (float)binsPerMag2 - minIndex2);

    if(valueInt1 < 0 && binsInLog1 == 1)
      valueInt1 = 0;
    if(valueInt2 < 0 && binsInLog2 == 1)
      valueInt2 = 0;
    
    assert(valueInt1 >= 0);
    assert(valueInt2 >= 0);
    assert(valueInt1 < numBins1);
    assert(valueInt2 < numBins2);
    
    histogramNum[valueInt1 * numBins2 + valueInt2] += numDensEndSnap[gal];
  }
  
  /* communicate histograms */
#ifdef MPI
  double *recvNumHistogram = allocate_array_double(numBins1 * numBins2, "recvNumHistogram");
  
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1*numBins2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
        histogramNum[valueInt1 * numBins2 + valueInt2] = recvNumHistogram[valueInt1 * numBins2 + valueInt2];
    }
  }
  
  free(recvNumHistogram);
#endif
  
  double binwidth1 = 0.;
  if(binsInLog1 == 1)
    binwidth1 = log10(values1[1]) - log10(values1[0]);
  else
    binwidth1 = values1[1] - values1[0];
  
  double binwidth2 = 0.;
  if(binsInLog2 == 1)
    binwidth2 = log10(values2[1]) - log10(values2[0]);
  else
    binwidth2 = values2[1] - values2[0];
  
  double sum = 0.;
  for(int i=0; i<numBins1; i++)
  {
    for(int j=0; j<numBins2; j++)
    {
      sum += histogramNum[i * numBins2 + j];
      histogramNum[i * numBins2 + j] = histogramNum[i * numBins2 + j] / (binwidth1 * binwidth2);
    }
  }
  
  /* write histograms */   
  write_2D_histogram_analytic(numBins1, values1, numBins2, values2, histogramNum, thisRank, filename);

  free(histogramNum);
  free(values1);
  free(values2);
}

void write_1D_histogram_analytic(int num_bins, float *values, double *hist, int thisRank, char *filename)
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
        for(int i=0; i<num_bins-1; i++)
        {
            printf("%e\t%e\t%e\n", values[i], values[i+1], hist[i]);
            fprintf(f, "%e\t%e\t%e\n", values[i], values[i+1], hist[i]);
        }
        
        fclose(f);
    }
}

void write_2D_histogram_analytic(int numBins1, float *values1, int numBins2, float *values2, double *numHistogram, int thisRank, char *filename)
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
    
    for(int valueInt1=0; valueInt1<numBins1; valueInt1++)
    {
      for(int valueInt2=0; valueInt2<numBins2; valueInt2++)
      {
        fprintf(f, "%e\t%e\t%e\n", values1[valueInt1], values2[valueInt2], numHistogram[valueInt1 * numBins2  + valueInt2]);
      }
    }
    
    fclose(f);
  }
}