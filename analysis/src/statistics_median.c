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
#include "find_median.h"
#include "statistics.h"
#include "statistics_median.h"

void calc_2D_histogram_history_median(int numGal, double *property, int numSnaps, float *redshifts, int currSnap, double *binProperty1, double *binProperty2, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, int size, char *filename)
{  
  int currSnap1 = 0, currSnap2 = 0;
  if(containsHistories1 == 1)
    currSnap1 = currSnap;
  if(containsHistories2 == 1)
    currSnap2 = currSnap;
  
  /* get minimum and maximum */
  float minValue1 = 0., maxValue1 = 0.;
  float minValue2 = 0., maxValue2 = 0.;
  float minValue = 0., maxValue = 0.;
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty1, &minValue1, &maxValue1, currSnap1);
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty2, &minValue2, &maxValue2, currSnap2);
  get_min_and_max_galaxy_property_with_histories(numGal, property, &minValue, &maxValue, currSnap);

  /* find binning */
  int numBins1 = 0, numBins2 = 0;
  float value1 = 0, value2 = 0;
  float *values1 = NULL, *values2 = NULL;
  int minIndex1 = 0, minIndex2 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  get_bins(minValue2, maxValue2, binsInLog2, binsPerMag2, &numBins2, &values2, &minIndex2);
  
  /* construct histogram */
  int valueInt1 = 0, valueInt2 = 0;
  int numAlloc = numGal / (numBins1 * numBins2);
  float *histogramMedian = allocate_array_float(numBins1 * numBins2 * (currSnap + 1), "histogram");
  float **histogram = allocate_array_float_pointer(numBins1 * numBins2 * (currSnap + 1), "histogram");
  int *histogramNum = allocate_array_int(numBins1 * numBins2 * (currSnap + 1), "histogramNum");
  int *histogramNumAlloc = allocate_array_int(numBins1 * numBins2 * (currSnap + 1), "histogramNumAlloc");
    
  for(int snap=0; snap<=currSnap; snap++)
  {
    for(valueInt1=0; valueInt1<numBins1; valueInt1++)
    {
      for(valueInt2=0; valueInt2<numBins2; valueInt2++)
      {
        histogramNumAlloc[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap] = numAlloc;
        histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap] = allocate_array_float(numAlloc, "histogram");
      }
    }
  }
  
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
        histogramNum[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap] += 1;
        if(histogramNum[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap] >= histogramNumAlloc[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap])
        {
          histogramNumAlloc[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap] += numAlloc;
          histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap] = realloc(histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap], histogramNumAlloc[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap]*sizeof(float));
        }
        histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap][histogramNum[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap]-1] = property[gal * (currSnap + 1) + snap];
      }
    }
  }
    
  for(int snap=0; snap<=currSnap; snap++)
  {
    for(valueInt1=0; valueInt1<numBins1; valueInt1++)
    {
      for(valueInt2=0; valueInt2<numBins2; valueInt2++)
      {
        /* cut down histogram to actual values */
        histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap] = realloc(histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap], histogramNum[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap]*sizeof(float));
        
        for(int i=0; i<histogramNum[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap]; i++)
        {
          if(histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap][i] < 0.)
            printf("Invalid value: %e\n", histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap][i]);
        }
        
        /* find median */
        histogramMedian[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap] = get_median(histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap], (long int)histogramNum[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap], thisRank, size, minValue, maxValue);
        
        free(histogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap]);
      }
    }
  }
  
  /* communicate histograms */
#ifdef MPI
  int *recvNumHistogram = allocate_array_int(numBins1 * numBins2 * (currSnap + 1), "recvNumHistogram");
  
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1*numBins2*(currSnap+1), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(int snap=0; snap<=currSnap; snap++)
  {
    for(valueInt1=0; valueInt1<numBins1; valueInt1++)
    {
      for(valueInt2=0; valueInt2<numBins2; valueInt2++)
      {
        histogramNum[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap] = recvNumHistogram[(valueInt1 * numBins2 + valueInt2) * (currSnap + 1) + snap];
      }
    }
  }
  
  free(recvNumHistogram);
#endif
          
  /* write histograms */        
  write_2D_histogram_history_median(numBins1, values1, numBins2, values2, currSnap, redshifts, histogramMedian, histogramNum, thisRank, filename);

  free(histogramMedian);
  free(histogram);
  free(histogramNum);
  free(histogramNumAlloc);
  free(values1);
  free(values2);
}

void calc_2D_histogram_median(int numGal, double *property, int currSnap, double *binProperty1, double *binProperty2, int binsInLog1, int binsInLog2, int binsPerMag1, int binsPerMag2, int containsHistories1, int containsHistories2, int thisRank, int size, char *filename)
{
  int currSnap1 = 0, currSnap2 = 0;
  if(containsHistories1 == 1)
    currSnap1 = currSnap;
  if(containsHistories2 == 1)
    currSnap2 = currSnap;
  
  /* get minimum and maximum */
  float minValue1 = 0., maxValue1 = 0.;
  float minValue2 = 0., maxValue2 = 0.;
  float minValue = 0., maxValue = 0.;
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty1, &minValue1, &maxValue1, currSnap1);
  get_min_and_max_galaxy_property_with_histories(numGal, binProperty2, &minValue2, &maxValue2, currSnap2);
  get_min_and_max_galaxy_property_with_histories(numGal, property, &minValue, &maxValue, 0);

  /* find binning */
  int numBins1 = 0, numBins2 = 0;
  float value1 = 0, value2 = 0;
  float *values1 = NULL, *values2 = NULL;
  int minIndex1 = 0, minIndex2 = 0;

  get_bins(minValue1, maxValue1, binsInLog1, binsPerMag1, &numBins1, &values1, &minIndex1);
  get_bins(minValue2, maxValue2, binsInLog2, binsPerMag2, &numBins2, &values2, &minIndex2);
  
  /* construct histogram */
  int valueInt1 = 0, valueInt2 = 0;
  int numAlloc = numGal / (numBins1 * numBins2);
  float *histogramMedian = allocate_array_float(numBins1 * numBins2, "histogram");
  float **histogram = allocate_array_float_pointer(numBins1 * numBins2, "histogram");
  int *histogramNum = allocate_array_int(numBins1 * numBins2, "histogramNum");
  int *histogramNumAlloc = allocate_array_int(numBins1 * numBins2, "histogramNumAlloc");
    
  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
      histogramNumAlloc[valueInt1 * numBins2 + valueInt2] = numAlloc;
      histogram[valueInt1 * numBins2 + valueInt2] = allocate_array_float(numAlloc, "histogram");
    }
  }
  
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
      histogramNum[valueInt1 * numBins2 + valueInt2] += 1;
      if(histogramNum[valueInt1 * numBins2 + valueInt2] >= histogramNumAlloc[valueInt1 * numBins2 + valueInt2])
      {
        histogramNumAlloc[valueInt1 * numBins2 + valueInt2] += numAlloc;
        histogram[valueInt1 * numBins2 + valueInt2] = realloc(histogram[valueInt1 * numBins2 + valueInt2], histogramNumAlloc[valueInt1 * numBins2 + valueInt2]*sizeof(float));
      }
      histogram[valueInt1 * numBins2 + valueInt2][histogramNum[valueInt1 * numBins2 + valueInt2]-1] = property[gal];
    }
  }
    
  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
      /* cut down histogram to actual values */
      histogram[valueInt1 * numBins2 + valueInt2] = realloc(histogram[valueInt1 * numBins2 + valueInt2], histogramNum[valueInt1 * numBins2 + valueInt2]*sizeof(float));
      
      /* find median */
      histogramMedian[valueInt1 * numBins2 + valueInt2] = get_median(histogram[valueInt1 * numBins2 + valueInt2], (long int)histogramNum[valueInt1 * numBins2 + valueInt2], thisRank, size, minValue, maxValue);
      
      free(histogram[valueInt1 * numBins2 + valueInt2]);
    }
  }
  
  /* communicate histograms */
#ifdef MPI
  int *recvNumHistogram = allocate_array_int(numBins1 * numBins2, "recvNumHistogram");
  
  MPI_Allreduce(histogramNum, recvNumHistogram, numBins1*numBins2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for(valueInt1=0; valueInt1<numBins1; valueInt1++)
  {
    for(valueInt2=0; valueInt2<numBins2; valueInt2++)
    {
      histogramNum[valueInt1 * numBins2 + valueInt2] = recvNumHistogram[valueInt1 * numBins2 + valueInt2];
    }
  }
  
  free(recvNumHistogram);
#endif
          
  /* write histograms */        
  write_2D_histogram_median(numBins1, values1, numBins2, values2, histogramMedian, histogramNum, thisRank, filename);

  free(histogramMedian);
  free(histogram);
  free(histogramNum);
  free(histogramNumAlloc);
  free(values1);
  free(values2);
}


void write_2D_histogram_history_median(int numBins1, float *values1, int numBins2, float *values2, int currSnap, float *redshifts, float *histogram, int *numHistogram, int thisRank, char *filename)
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

void write_2D_histogram_median(int numBins1, float *values1, int numBins2, float *values2, float *histogram, int *numHistogram, int thisRank, char *filename)
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
