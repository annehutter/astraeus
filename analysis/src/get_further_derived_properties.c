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
#include "dconfObj.h"
#include "outgal.h"
#include "derive_properties.h"
#include "operations_on_properties.h"
#include "ion_emissivity.h"
#include "SNfeedback.h"
#include "get_further_derived_properties.h"

void prepare_stellarmasshistory(double *stellarmasshistory, int numGal, int currSnap)
{
  for(int gal=0; gal<numGal; gal++)
  {
    for(int snap=0; snap<=currSnap; snap++)
    {
      if(stellarmasshistory[gal*(currSnap + 1) + snap] < 0.)
        stellarmasshistory[gal*(currSnap + 1) + snap] = 0.;
    }
  }
}

int check_derived_property_present(char *propertyName)
{
  int result = 0;
  
  /* check whether a derived property is being asked for */
  if(strcmp(propertyName, "fej") == 0 || strcmp(propertyName, "fesc") == 0 || strcmp(propertyName, "fescFej") == 0 || strcmp(propertyName, "Nion") == 0)
    result = 1;
  
  return result;
}

double get_max_of_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  int numGal = 0;
  double *thisProperty = getThisProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, propertyName, currSnap, times, &numGal);

  assert(numGal == get_numGal_at_snap(sizeListEquals, numTrees));
  
  double max = 0.;
  
  for(int gal=0; gal<numGal; gal++)
  {
    if(thisProperty[gal] > max)
      max = thisProperty[gal];
  }
  
  free(thisProperty);
  
  return max;
}

double *get_derived_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  double *thisProperty = NULL;
  
  if(strcmp(propertyName,"fej") == 0)
  {
    thisProperty = get_SNejection_property(simParam, treeList, numTrees, listEquals, index, sizeListEquals, currSnap, times);
  }
  else if(strcmp(propertyName,"fescFej") == 0)
  {
    thisProperty = get_fesc_SNejection_property(simParam, treeList, numTrees, listEquals, index, sizeListEquals, currSnap, times);
  }
  else if(strcmp(propertyName,"fesc") == 0)
  {
    thisProperty = get_fesc_property(simParam, numTrees, sizeListEquals);
  }
  else if(strcmp(propertyName, "Nion") == 0)
  {
    thisProperty = get_ionEmissivity_property(simParam, treeList, numTrees, listEquals, index, sizeListEquals, currSnap, times);
  }
  
  return thisProperty;
}

double *get_ionEmissivity_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times)
{
  int numGal = 0; 
  
  double *thisStellarMassHistory = getThisPropertyHistory(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "Mstar", currSnap, times, &numGal);
  assert(numGal == get_numGal_at_snap(sizeListEquals, numTrees));

  double *Mvir = getThisProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "Mvir", currSnap, times, &numGal);
  assert(numGal == get_numGal_at_snap(sizeListEquals, numTrees));

  prepare_stellarmasshistory(thisStellarMassHistory, numGal, currSnap);
  
  double *thisProperty = allocate_array_double(numGal, "ionEmissivity");
  
  for(int gal=0; gal<numGal; gal++)
  {
    thisProperty[gal] = (double)(1.e-50*get_nion_sps(&(thisStellarMassHistory[gal*(currSnap+1)]), currSnap, simParam));
  }
  
  free(Mvir);
  free(thisStellarMassHistory);
  
  return thisProperty;
}

double *get_SNejection_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times)
{
  int numGal = 0;
  double *thisProperty = NULL;
  double *thisStellarMassHistory = NULL;
  double *thisMgasIni = NULL;
  double *thisMvirDivRvir = NULL;

  /* get star formation history */
  thisStellarMassHistory = getThisPropertyHistory(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "Mstar", currSnap, times, &numGal);
  assert(numGal == get_numGal_at_snap(sizeListEquals, numTrees));
  prepare_stellarmasshistory(thisStellarMassHistory, numGal, currSnap);

  thisMgasIni = getThisProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "MgasIni", currSnap, times, &numGal);
  assert(numGal == get_numGal_at_snap(sizeListEquals, numTrees));

  thisMvirDivRvir = getThisProperty(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "Mvir%Rvir", currSnap, times, &numGal);
  assert(numGal == get_numGal_at_snap(sizeListEquals, numTrees));

  /* get ejection efficiency */
  thisProperty = get_SNejection_efficiency(simParam, numGal, thisMvirDivRvir, thisMgasIni, thisStellarMassHistory, currSnap);

  free(thisStellarMassHistory);
  free(thisMgasIni);
  free(thisMvirDivRvir);
  
  return thisProperty;
}

double *get_fesc_SNejection_property(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, float *times)
{
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
  double *SNejectionEfficiency = get_SNejection_property(simParam, treeList, numTrees, listEquals, index, sizeListEquals, currSnap, times);
  
  double feff = get_max_of_property(simParam, treeList, numTrees, listEquals, index, sizeListEquals, "feff", currSnap, times);
  
  double *fesc = allocate_array_double(numGal, "fesc");
  
  for(int gal=0; gal<numGal; gal++)
  {
    fesc[gal] = 1.;
    if(SNejectionEfficiency[gal] > 0.)
      fesc[gal] = feff/SNejectionEfficiency[gal];
    if(fesc[gal] > 1.)
      fesc[gal] = 1.;
    fesc[gal] = fesc[gal] * simParam->fesc;
  }
  
  free(SNejectionEfficiency);
  
  return fesc;
}

double *get_fesc_property(dconfObj_t simParam, int numTrees, int *sizeListEquals)
{
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
  double *fesc = allocate_array_double(numGal, "fesc");
    
  for(int gal=0; gal<numGal; gal++)
  {
    fesc[gal] = simParam->fesc;
  }
    
  return fesc;
}