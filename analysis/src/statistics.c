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

int get_max(int *arrayInt, int length)
{
  int max = -10;
  for(int i=0; i<length; i++)
    if(arrayInt[i] > max)
      max = arrayInt[i];
  
  return max;
}

void get_min_and_max_galaxy_property(int numGal, double *property, float *min, float *max)
{
  float minValue = 1.e60;
  float maxValue = -1.e60;
  
  for(int gal=0; gal<numGal; gal++)
  {
    if(property[gal] < minValue)
      minValue = property[gal];
    if(property[gal] > maxValue)
      maxValue = property[gal];
  }
  
#ifdef MPI
  float recvMinValue = 0.;
  float recvMaxValue = 0.;
    
  MPI_Allreduce(&minValue, &recvMinValue, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&maxValue, &recvMaxValue, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  
  minValue = recvMinValue;
  maxValue = recvMaxValue;
#endif
  
  *min = minValue;
  *max = maxValue;
}

void get_min_and_max_galaxy_property_with_histories(int numGal, double *property, float *min, float *max, int currSnap)
{
  float minValue = 1.e60;
  float maxValue = -1.e60;
  
  for(int gal=0; gal<numGal; gal++)
  {
    if(property[gal * (currSnap + 1) + currSnap] < minValue)
      minValue = property[gal * (currSnap + 1) + currSnap];
    if(property[gal * (currSnap + 1) + currSnap] > maxValue)
      maxValue = property[gal * (currSnap + 1) + currSnap];
  }
  
#ifdef MPI
  float recvMinValue = 0.;
  float recvMaxValue = 0.;
    
  MPI_Allreduce(&minValue, &recvMinValue, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&maxValue, &recvMaxValue, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
  
  minValue = recvMinValue;
  maxValue = recvMaxValue;
#endif
  
  *min = minValue;
  *max = maxValue;
}

void get_bins(float minValue, float maxValue, int binsInLog, int binsPerMag, int *numBinsIn, float **valuesIn, int *minIndexIn)
{
  int numBins = 0;
  int minIndex = 0;
  float *values = NULL;
  
  if(binsInLog == 1)
  {
    if(minValue <= 0.) minValue = maxValue * 1.e-10; 
    numBins = (int)( (log10(maxValue) - log10(minValue)) * (float)binsPerMag ) + 2;
    minIndex = floor(log10(minValue) * (float)binsPerMag);
    
    if(numBins > 0) values = allocate_array_float(numBins, "values");
    for(int i=0; i<numBins; i++)
    {
      values[i] = pow(10., (minIndex + (float)i) / (float)binsPerMag);
    }
  }
  else
  {
    numBins = (int)( (maxValue - minValue) * (float)binsPerMag ) + 2;
    minIndex = floor(minValue * (float)binsPerMag);
    
    if(numBins > 0) values = allocate_array_float(numBins, "values");
    for(int i=0; i<numBins; i++)
    {
      values[i] = (minIndex + (float)i) / (float)binsPerMag;
    }
  }
  
  *numBinsIn = numBins;
  *valuesIn = values;
  *minIndexIn = minIndex;
}

void calc_1D_histogram_history(int numGal, double *property, int numSnaps, float *redshifts, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1, int thisRank, char *filename)
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
  float *histogram = allocate_array_float(numBins1 * (currSnap + 1), "histogram");
  int *histogramNum = allocate_array_int(numBins1 * (currSnap + 1), "histogramNum");
  
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
    
    for(int snap=0; snap<=currSnap; snap++)
    {
      if(property[gal * (currSnap + 1) + snap] >= 0.)
      {
        histogram[valueInt1 * (currSnap + 1) + snap] += property[gal * (currSnap + 1) + snap];
        histogramNum[valueInt1 * (currSnap + 1) + snap] += 1;
      }
    }
  }
  
  /* communicate histograms */
#ifdef MPI
  float *recvHistogram = allocate_array_float(numBins1 * (currSnap + 1), "recvHistogram");
  int *recvNumHistogram = allocate_array_int(numBins1 * (currSnap + 1), "recvNumHistogram");
  
  MPI_Allreduce(histogram, recvHistogram, numBins1*(currSnap + 1), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1*(currSnap + 1), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(int snap=0; snap<=currSnap; snap++)
    {
      histogram[valueInt1 * (currSnap + 1) + snap] = recvHistogram[valueInt1 * (currSnap + 1) + snap];
      histogramNum[valueInt1 * (currSnap + 1) + snap] = recvNumHistogram[valueInt1 * (currSnap + 1) + snap];
    }
  }
  
  free(recvHistogram);
  free(recvNumHistogram);
#endif

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(int snap=0; snap<=currSnap; snap++)
    {
        if(histogramNum[valueInt1 * (currSnap + 1) + snap] > 0)
          histogram[valueInt1 * (currSnap + 1) + snap] /= (float) histogramNum[valueInt1 * (currSnap + 1) + snap];
    }
  }
          
  /* write histograms */        
  write_1D_histogram_history(numBins1, values1, currSnap, redshifts, histogram, histogramNum, thisRank, filename);

  free(histogram);
  free(histogramNum);
  free(values1);
}

void calc_2D_histogram_history(int numGal, double *property, int numSnaps, float *redshifts, int currSnap, double *binProperty1, double *binProperty2, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, char *filename)
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
  
  /* construct histogram */
  int valueInt1 = 0, valueInt2 = 0;
  float *histogram = allocate_array_float(numBins1 * numBins2 * (currSnap + 1), "histogram");
  int *histogramNum = allocate_array_int(numBins1 * numBins2 * (currSnap + 1), "histogramNum");
   
  for(int gal=0; gal<numGal; gal++)
  {
    if(binsInLog1 == 1)
      value1 = log10(binProperty1[gal * (currSnap1 + 1) + currSnap1]);
    else
      value1 = binProperty1[gal * (currSnap1 + 1) + currSnap1];
    
    if(binsInLog2 == 1)
      value2 = log10(binProperty2[gal * (currSnap2 + 1) + currSnap2]);
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
    
    for(int snap=0; snap<=currSnap; snap++)
    {
      if(property[gal * (currSnap + 1) + snap] >= 0.)
      {
        histogram[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap] += property[gal * (currSnap + 1) + snap];
        histogramNum[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap] += 1;
      }
    }
  }
  
  /* communicate histograms */
#ifdef MPI
  float *recvHistogram = allocate_array_float(numBins1 * numBins2 * (currSnap + 1), "recvHistogram");
  int *recvNumHistogram = allocate_array_int(numBins1 * numBins2 * (currSnap + 1), "recvNumHistogram");
  
  MPI_Allreduce(histogram, recvHistogram, numBins1*numBins2*(currSnap + 1), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1*numBins2*(currSnap + 1), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
      for(int snap=0; snap<=currSnap; snap++)
      {
        histogram[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap] = recvHistogram[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap];
        histogramNum[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap] = recvNumHistogram[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap];
      }
    }
  }
  
  free(recvHistogram);
  free(recvNumHistogram);
#endif

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
      for(int snap=0; snap<=currSnap; snap++)
      {
          if(histogramNum[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap] > 0)
            histogram[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap] /= (float) histogramNum[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap];
      }
    }
  }
          
  /* write histograms */        
  write_2D_histogram_history(numBins1, values1, numBins2, values2, currSnap, redshifts, histogram, histogramNum, thisRank, filename);

  free(histogram);
  free(histogramNum);
  free(values1);
  free(values2);
}

void write_1D_histogram_history(int numBins1, float *values1, int currSnap, float *redshifts, float *histogram, int *numHistogram, int thisRank, char *filename)
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
      for(int snap=0; snap<=currSnap; snap++)
      {
        fprintf(f, "%e\t%e\t%e\t%d\n", values1[valueInt1], redshifts[snap], histogram[valueInt1 * (currSnap + 1) + snap], numHistogram[valueInt1 * (currSnap + 1) + snap]);
      }
    }
    
    fclose(f);
  }
}

void write_2D_histogram_history(int numBins1, float *values1, int numBins2, float *values2, int currSnap, float *redshifts, float *histogram, int *numHistogram, int thisRank, char *filename)
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
        for(int snap=0; snap<=currSnap; snap++)
        {
          fprintf(f, "%e\t%e\t%e\t%e\t%d\n", values1[valueInt1], values2[valueInt2], redshifts[snap], histogram[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap], numHistogram[valueInt1 * numBins2 * (currSnap + 1) + valueInt2 * (currSnap + 1) + snap]);
        }
      }
    }
    
    fclose(f);
  }
}

void calc_1D_histogram_value(int numGal, double *property, int currSnap, double *binProperty1, int binsInLog1, int binsPerMag1, int containsHistories1, int thisRank, char *filename)
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
  float *histogram = allocate_array_float(numBins1, "histogram");
  int *histogramNum = allocate_array_int(numBins1, "histogramNum");
  
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
    
    if(property[gal] >= 0.)
    {
      histogram[valueInt1] += property[gal];
      histogramNum[valueInt1] += 1;
    }
  }
  
  /* communicate histograms */
#ifdef MPI
  float *recvHistogram = allocate_array_float(numBins1, "recvHistogram");
  int *recvNumHistogram = allocate_array_int(numBins1, "recvNumHistogram");
  
  MPI_Allreduce(histogram, recvHistogram, numBins1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    histogram[valueInt1] = recvHistogram[valueInt1];
    histogramNum[valueInt1] = recvNumHistogram[valueInt1];
  }
  
  free(recvHistogram);
  free(recvNumHistogram);
#endif

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    if(histogramNum[valueInt1] > 0)
      histogram[valueInt1] /= (float) histogramNum[valueInt1];
  }
          
  /* write histograms */        
  write_1D_histogram_value(numBins1, values1, histogram, histogramNum, thisRank, filename);

  free(histogram);
  free(histogramNum);
  free(values1);
}

void calc_2D_histogram_value(int numGal, double *property, int currSnap, double *binProperty1, double *binProperty2, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, char *filename)
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
    
  /* construct histogram */
  int valueInt1 = 0, valueInt2 = 0;
  float *histogram = allocate_array_float(numBins1 * numBins2, "histogram");
  int *histogramNum = allocate_array_int(numBins1 * numBins2, "histogramNum");
   
  for(int gal=0; gal<numGal; gal++)
  {
    if(binsInLog1 == 1)
      value1 = log10(binProperty1[gal * (currSnap1 + 1) + currSnap1]);
    else
      value1 = binProperty1[gal * (currSnap1 + 1) + currSnap1];
    
    if(binsInLog2 == 1)
      value2 = log10(binProperty2[gal * (currSnap2 + 1) + currSnap2]);
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
    
    if(property[gal] >= 0.)
    {
      histogram[valueInt1 * numBins2 + valueInt2] += property[gal];
      histogramNum[valueInt1 * numBins2 + valueInt2] += 1;
    }
  }
  
  /* communicate histograms */
#ifdef MPI
  float *recvHistogram = allocate_array_float(numBins1 * numBins2, "recvHistogram");
  int *recvNumHistogram = allocate_array_int(numBins1 * numBins2, "recvNumHistogram");
  
  MPI_Allreduce(histogram, recvHistogram, numBins1*numBins2, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1*numBins2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
      histogram[valueInt1 * numBins2 + valueInt2] = recvHistogram[valueInt1 * numBins2 + valueInt2];
      histogramNum[valueInt1 * numBins2 + valueInt2] = recvNumHistogram[valueInt1 * numBins2 + valueInt2];
    }
  }
  
  free(recvHistogram);
  free(recvNumHistogram);
#endif

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
      if(histogramNum[valueInt1 * numBins2 + valueInt2] > 0)
        histogram[valueInt1 * numBins2 + valueInt2] /= (float) histogramNum[valueInt1 * numBins2 + valueInt2];
    }
  }
          
  /* write histograms */        
  write_2D_histogram_value(numBins1, values1, numBins2, values2, histogram, histogramNum, thisRank, filename);

  free(histogram);
  free(histogramNum);
  free(values1);
  free(values2);
}

void calc_3D_histogram_value(int numGal, double *property, int currSnap, double *binProperty1, double *binProperty2, double *binProperty3, int binsInLog1, int binsInLog2, int binsInLog3, int binsPerMag1, int binsPerMag2, int binsPerMag3, int containsHistories1, int containsHistories2, int containsHistories3, int thisRank, char *filename)
{
  int currSnap1 = 0, currSnap2 = 0, currSnap3 = 0;
  if(containsHistories1 == 1)
    currSnap1 = currSnap;
  if(containsHistories2 == 1)
    currSnap2 = currSnap;
  if(containsHistories3 == 1)
    currSnap3 = currSnap;
  
  /* get minimum and maximum */
  float minValue1 = 0., maxValue1 = 0.;
  float minValue2 = 0., maxValue2 = 0.;
  float minValue3 = 0., maxValue3 = 0.;
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty1, &minValue1, &maxValue1, currSnap1);
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty2, &minValue2, &maxValue2, currSnap2);
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty3, &minValue3, &maxValue3, currSnap3);

  /* find binning */
  int numBins1 = 0, numBins2 = 0, numBins3 = 0;
  float value1 = 0, value2 = 0, value3 = 0;
  float *values1 = NULL, *values2 = NULL, *values3 = NULL;
  int minIndex1 = 0, minIndex2 = 0, minIndex3 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  get_bins(minValue2, maxValue2, binsInLog2, binsPerMag2, &numBins2, &values2, &minIndex2);
  get_bins(minValue3, maxValue3, binsInLog3, binsPerMag3, &numBins3, &values3, &minIndex3);
  
  /* construct histogram */
  int valueInt1 = 0, valueInt2 = 0, valueInt3 = 0;
  float *histogram = allocate_array_float(numBins1 * numBins2 * numBins3, "histogram");
  int *histogramNum = allocate_array_int(numBins1 * numBins2 * numBins3, "histogramNum");
   
  for(int gal=0; gal<numGal; gal++)
  {
    if(binsInLog1 == 1)
      value1 = log10(binProperty1[gal * (currSnap1 + 1) + currSnap1]);
    else
      value1 = binProperty1[gal * (currSnap1 + 1) + currSnap1];
    
    if(binsInLog2 == 1)
      value2 = log10(binProperty2[gal * (currSnap2 + 1) + currSnap2]);
    else
      value2 = binProperty2[gal * (currSnap2 + 1) + currSnap2];
        
    if(binsInLog3 == 1)
      value3 = log10(binProperty3[gal * (currSnap3 + 1) + currSnap3]);
    else
      value3 = binProperty3[gal * (currSnap3 + 1) + currSnap3];
    
    valueInt1 = floor(value1 * (float)binsPerMag1 - minIndex1);
    valueInt2 = floor(value2 * (float)binsPerMag2 - minIndex2);
    valueInt3 = floor(value3 * (float)binsPerMag3 - minIndex3);

    if(valueInt1 < 0 && binsInLog1 == 1)
      valueInt1 = 0;
    if(valueInt2 < 0 && binsInLog2 == 1)
      valueInt2 = 0;
    if(valueInt3 < 0 && binsInLog3 == 1)
      valueInt3 = 0;
    
    assert(valueInt1 >= 0);
    assert(valueInt2 >= 0);
    assert(valueInt3 >= 0);
    assert(valueInt1 < numBins1);
    assert(valueInt2 < numBins2);
    assert(valueInt3 < numBins3);
    
    if(property[gal] >= 0.)
    {
      histogram[valueInt1 * numBins2 * numBins3 + valueInt2 * numBins3 + valueInt3] += property[gal];
      histogramNum[valueInt1 * numBins2 *numBins3 + valueInt2 * numBins3 + valueInt3] += 1;
    }
  }
  
  /* communicate histograms */
#ifdef MPI
  float *recvHistogram = allocate_array_float(numBins1 * numBins2 * numBins3, "recvHistogram");
  int *recvNumHistogram = allocate_array_int(numBins1 * numBins2 * numBins3, "recvNumHistogram");
  
  MPI_Allreduce(histogram, recvHistogram, numBins1*numBins2*numBins3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1*numBins2*numBins3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
      for(valueInt3=0; valueInt3<numBins3; valueInt3++)
      {
        histogram[valueInt1 * numBins2 *numBins3 + valueInt2 * numBins3 + valueInt3] = recvHistogram[valueInt1 * numBins2 * numBins3 + valueInt2 * numBins3 + valueInt3];
        histogramNum[valueInt1 * numBins2 * numBins3 + valueInt2 * numBins3 + valueInt3] = recvNumHistogram[valueInt1 * numBins2 * numBins3 + valueInt2 * numBins3 + valueInt3];
      }
    }
  }
  
  free(recvHistogram);
  free(recvNumHistogram);
#endif

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
      for(valueInt3=0; valueInt3<numBins3; valueInt3++)
      {
        if(histogramNum[valueInt1 * numBins2 * numBins3 + valueInt2 * numBins3 + valueInt3] > 0)
          histogram[valueInt1 * numBins2 * numBins3 + valueInt2 * numBins3 + valueInt3] /= (float) histogramNum[valueInt1 * numBins2 * numBins3 + valueInt2 * numBins3 + valueInt3];
      }
    }
  }
          
  /* write histograms */        
  write_3D_histogram_value(numBins1, values1, numBins2, values2, numBins3, values3, histogram, histogramNum, thisRank, filename);

  free(histogram);
  free(histogramNum);
  free(values1);
  free(values2);
  free(values3);
}

void write_1D_histogram_value(int numBins1, float *values1, float *histogram, int *numHistogram, int thisRank, char *filename)
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
      fprintf(f, "%e\t%e\t%d\n", values1[valueInt1], histogram[valueInt1], numHistogram[valueInt1]);
    }
    
    fclose(f);
  }
}

void write_2D_histogram_value(int numBins1, float *values1, int numBins2, float *values2, float *histogram, int *numHistogram, int thisRank, char *filename)
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
        fprintf(f, "%e\t%e\t%e\t%d\n", values1[valueInt1], values2[valueInt2], histogram[valueInt1 * numBins2 + valueInt2], numHistogram[valueInt1 * numBins2 + valueInt2]);
      }
    }
    
    fclose(f);
  }
}

void write_3D_histogram_value(int numBins1, float *values1, int numBins2, float *values2, int numBins3, float *values3, float *histogram, int *numHistogram, int thisRank, char *filename)
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
        for(int valueInt3=0; valueInt3<numBins3; valueInt3++)
        {
          fprintf(f, "%e\t%e\t%e\t%e\t%d\n", values1[valueInt1], values2[valueInt2], values3[valueInt3], histogram[valueInt1 * numBins2 * numBins3 + valueInt2 * numBins3 + valueInt3], numHistogram[valueInt1 * numBins2 * numBins3 + valueInt2 * numBins3 + valueInt3]);
        }
      }
    }
    
    fclose(f);
  }
}

void calc_1D_histogram(int numGal, double *binProperty1, int currSnap, int binsInLog1, int binsPerMag1, int containsHistories1, float volumeInMpc, int cumulative, int thisRank, char *filename)
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
      
    histogramNum[valueInt1] += 1.;
  }
  
  double sumNum = 0;
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
  if(thisRank == 0)
    printf("1D histogram: numGal = %d\t numSum = %d\n", numGal, (int)sumNum);
  
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
      sum += histogramNum[i] / (volumeInMpc * binwidth);
      histogramNum[i] = sum;
    }
  }
  else{
    for(int i=0; i<numBins1; i++)
    {
      histogramNum[i] = histogramNum[i] / (volumeInMpc * binwidth);
    }
  }
      
  write_1D_histogram(numBins1, values1, histogramNum, thisRank, filename);
      
  free(values1);
  free(histogramNum);
}

void calc_2D_histogram(int numGal, double *binProperty1, double *binProperty2, int currSnap, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, float volumeInMpc, int thisRank, char *filename)
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
    printf("2D histogram: minValue1 = %e\t maxValue1 = %e\t numBins1 = %d\n", minValue1, maxValue1, numBins1);
    printf("2D histogram: minValue2 = %e\t maxValue2 = %e\t numBins2 = %d\n", minValue2, maxValue2, numBins2);
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
    
    histogramNum[valueInt1 * numBins2 + valueInt2] += 1.;
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
      histogramNum[i * numBins2 + j] = histogramNum[i * numBins2 + j] / (volumeInMpc * binwidth1 * binwidth2);
    }
  }
  
  if(thisRank == 0)
    printf("2D histogram: numGal = %d\n", (int)sum);
  
  /* write histograms */   
  write_2D_histogram(numBins1, values1, numBins2, values2, histogramNum, thisRank, filename);

  free(histogramNum);
  free(values1);
  free(values2);
}

void write_1D_histogram(int num_bins, float *values, double *hist, int thisRank, char *filename)
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

void write_2D_histogram(int numBins1, float *values1, int numBins2, float *values2, double *numHistogram, int thisRank, char *filename)
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