#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#ifdef MPI
#include <mpi.h>
#endif

#include "const.h"
#include "utils.h"
#include "dconfObj.h"
#include "gal_gtree.h"
#include "metallicity.h"

#if defined WITHMETALS

metal_t *initMetal(int32_t endSnap)
{
  metal_t *thisMetal = malloc(sizeof(metal_t));
  
  if(thisMetal == NULL)
  {
    fprintf(stderr, "Could not allocate metal_t\n");
    exit(EXIT_FAILURE);
  }
  
  thisMetal->endSnap = endSnap;
  thisMetal->meanMetallicityIgm = allocate_array_double(endSnap, "meanMetallicityIgm");
  thisMetal->meanMetallicityIgmO = allocate_array_double(endSnap, "meanMetallicityIgmO");
  thisMetal->meanMetallicityIgmFe = allocate_array_double(endSnap, "meanMetallicityIgmFe");
  thisMetal->meanDustfractionIgm = allocate_array_double(endSnap, "meanDustfractionIgm");

  return thisMetal;
}

void deallocate_metal(metal_t *thisMetal)
{
  if(thisMetal->meanMetallicityIgm != NULL) free(thisMetal->meanMetallicityIgm);
  if(thisMetal->meanMetallicityIgmO != NULL) free(thisMetal->meanMetallicityIgmO);
  if(thisMetal->meanMetallicityIgmFe != NULL) free(thisMetal->meanMetallicityIgmFe);
  if(thisMetal->meanDustfractionIgm != NULL) free(thisMetal->meanDustfractionIgm);

  free(thisMetal);
}


int read_array_from_txtfile(char *fileName, double **thisArray)
{
  double *array = allocate_array_double(MAXLENGTH, "array");
  int numEntries = 0;
  int numAllocated = MAXLENGTH;
  
  FILE *f = NULL;
  
  if( (f=fopen(fileName, "r")) == NULL)
  {
    fprintf(stderr, "Could not open file\n");
    exit(EXIT_FAILURE);
  }
  
  while(fscanf(f, "%lf", &(array[numEntries])) == 1)
  {
    numEntries++;
    if(numEntries >= numAllocated)
    {
      numAllocated = numAllocated * 2;
      array = realloc(array, sizeof(double) * numAllocated);
    }
  }
  array = realloc(array, sizeof(double) * numEntries);
  
  fclose(f);
  
  *thisArray = array;
  
  return numEntries;
}

int where_loop(double array[], int n, double value)
{
  int index = 0;
  double diff, temp_diff;
  
  if(value <= array[0]) 
    return 0;
  else if(value >= array[n-1]) 
    return (n-1);
  else 
  {
    temp_diff = fabs(value - array[0]);
    index = 0;
    for(int i=1; i<n; i++)
    {
      diff = fabs(value - array[i]);
      if(diff <= temp_diff) 
      {
        temp_diff = diff;
        index = i;
      }
    }
    return index;
  }
}

/* This function is based on binary research tree */
int where(double array[], int n, double search)
{
  int first  = 0;
  int last   = n - 1;
  int middle = (first + last) / 2;
  
  while (first <= last)
  {
    if (array[middle] < search)
      first = middle + 1;
    else if (array[middle] == search)
    {
      return middle;
      break;
    }
    else
      last = middle - 1;

    middle = (first + last)/2;
  }
  if (first > last) 
  {
    if( fabs(search-array[middle]) > fabs(search-array[middle+1]) ) { middle = middle + 1;}
  }
  
  return middle;
}

/* Linear interpolation useful to compute the sfrAtFormation array */
double linear_interp(double x[], double y[], int length, double xx)
{
  double yy = 0;
  for(int i=0; i<(length-1); i++)
  {
    if((xx >= x[i]) & (xx <= x[i+1]))
    {
      yy = y[i] + ((y[i+1] - y[i]) / (x[i+1] - x[i])) * (xx - x[i]);
      break;
    }
  }
  return yy;
}

void add_ejected_metals_to_igm(gal_t *thisGal, int snap, double invTotMgas, metal_t *thisMetal)
{
  thisMetal->meanMetallicityIgm[snap] += thisGal->MmetalEj[0] * invTotMgas;
  thisMetal->meanMetallicityIgmO[snap] += thisGal->MmetalEj[1] * invTotMgas;
  thisMetal->meanMetallicityIgmFe[snap] += thisGal->MmetalEj[2] * invTotMgas;
  thisMetal->meanDustfractionIgm[snap] += thisGal->MdustEj * invTotMgas;
}

void add_igm_metals_to_gal(metal_t *thisMetal, gal_t *thisGal, int snap)
{
  thisGal->igmMetallicity[0] = thisMetal->meanMetallicityIgm[snap];
  thisGal->igmMetallicity[1] = thisMetal->meanMetallicityIgmO[snap];
  thisGal->igmMetallicity[2] = thisMetal->meanMetallicityIgmFe[snap];
  thisGal->igmDustFraction = thisMetal->meanDustfractionIgm[snap];
}

/* function to get the position of galaxies in next snapshot */
void map_metallicityIgm_to_gal(int numGtrees, gtree_t ***thisGtreeList, int minSnap, int maxSnap, metal_t *thisMetal)
{
  gtree_t **theseGtrees = *thisGtreeList;
  gtree_t *thisGtree = NULL;
  gal_t *thisGal;
  int32_t gal = 0;
  int32_t snapGal = 0;
  int32_t numGal = 0;
     
  /* LOOP OVER TREES */
  for(int gtree=0; gtree<numGtrees; gtree++)
  {
    thisGtree = theseGtrees[gtree];
    numGal = thisGtree->numGal;
    gal = thisGtree->walker;
    
    if(gal < numGal) 
    {
      snapGal = thisGtree->galaxies[gal].snapnumber;

      /* LOOP OVER SNAPSHOTS */
      for(int snap=minSnap+1; snap<maxSnap+1; snap++)
      {
        if(snapGal == snap)
        {
          /* LOOP OVER GALAXIES */
          thisGal = &(thisGtree->galaxies[gal]);
          while(thisGal->snapnumber == snap)
          { 
            add_igm_metals_to_gal(thisMetal, thisGal, snap-1);

            gal++;
            
            /* stop when reaching the end of the tree */
            if(gal == numGal)
              break;
            
            /* go to the next galaxy */
            thisGal = &(thisGtree->galaxies[gal]);
          }
          if(gal < numGal)
            snapGal = thisGtree->galaxies[gal].snapnumber;
        }
      }
    }
  }
}

/* Useful to update the metallicityIgm value (average */
void distribute_metallicityIgm(metal_t *thisMetal, int snap)
{
  if(snap != (thisMetal->endSnap-1))
  {
    thisMetal->meanMetallicityIgm[snap+1]   = thisMetal->meanMetallicityIgm[snap];
    thisMetal->meanMetallicityIgmO[snap+1]  = thisMetal->meanMetallicityIgmO[snap];
    thisMetal->meanMetallicityIgmFe[snap+1] = thisMetal->meanMetallicityIgmFe[snap]; 
    thisMetal->meanDustfractionIgm[snap+1]  = thisMetal->meanDustfractionIgm[snap];
  }
  
#ifdef MPI
  double recvMetallicityIgm, recvMetallicityIgmO, recvMetallicityIgmFe, recvDustfractionIgm;
  MPI_Allreduce(&(thisMetal->meanMetallicityIgm[snap]),  &recvMetallicityIgm,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  thisMetal->meanMetallicityIgm[snap] = recvMetallicityIgm;
  MPI_Allreduce(&(thisMetal->meanMetallicityIgmO[snap]), &recvMetallicityIgmO, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  thisMetal->meanMetallicityIgmO[snap] = recvMetallicityIgmO;
  MPI_Allreduce(&(thisMetal->meanMetallicityIgmFe[snap]), &recvMetallicityIgmFe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  thisMetal->meanMetallicityIgmFe[snap] = recvMetallicityIgmFe;
  MPI_Allreduce(&(thisMetal->meanDustfractionIgm[snap]), &recvDustfractionIgm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  thisMetal->meanDustfractionIgm[snap] = recvDustfractionIgm;
#endif
}

int *compute_indexes(int lower, int upper, int numMetallicityBins, double *metallicity, double *metallicityAtFormation)
{
  int *indexes = (int*) malloc((upper-lower) * sizeof(int));
  int j = 0;
  // Computing the integrals with the trapezoidal rule
  for(int i=(lower+1); i<upper; i++)
  {
    indexes[j] = where(metallicity, numMetallicityBins, metallicityAtFormation[i]);
    j += 1;
  }

  return indexes;
}

double compute_yield(double *yields, double *numStarsAtFormation, int numMassBins, int *indexes, int index_low, int index_up, int lowerIntegrationLimit, int upperIntegrationLimit, double stepsizeIntegration)
{
  double integral = yields[index_low * numMassBins + lowerIntegrationLimit] * numStarsAtFormation[lowerIntegrationLimit];
  int j = 0;
  // Computing the integrals with the trapezoidal rule
  for(int i=(lowerIntegrationLimit+1); i<upperIntegrationLimit; i++)
  {
    integral += 2 * (yields[indexes[j] * numMassBins + i] * numStarsAtFormation[i]);
    j += 1;
  }
  // index_up is already computed above
  integral += yields[index_up * numMassBins + upperIntegrationLimit] * numStarsAtFormation[upperIntegrationLimit];
  // 0.5 * stepsizeIntegration is dx/2 in the following
  integral  = 0.5 * stepsizeIntegration * integral;
    
  return integral;
}

double compute(double yield, double *numStarsAtFormation, int lowerIntegrationLimit, int upperIntegrationLimit, double stepsizeIntegration)
{
  double integral = numStarsAtFormation[lowerIntegrationLimit];
  for(int i=(lowerIntegrationLimit+1); i<upperIntegrationLimit; i++) 
  {
    integral += 2 * numStarsAtFormation[i];
  }
  integral += numStarsAtFormation[upperIntegrationLimit];
  // 0.5 * stepsizeIntegration is dx/2 in the following
  integral = 0.5 * stepsizeIntegration * yield * integral;
  
  return integral;
}

void compute_metal_mass(gal_t *thisGal, dconfObj_t simParam, float sfr, float deltat, float *G_t, float *e_Z, float *e_Z_O, float *e_Z_Fe, float *y_d) 
{
  int numMassBins = simParam->metalTables_numMassBins;
  int numMetallicityBins = simParam->metalTables_numMetallicityBins;
  double time = 0.;                     /* time corresponding to the current snapshot */
  double sfrAtFormation[numMassBins];            /* useful to define the array psi(t - lifetime) in the GCE equation */
  double metallicityAtFormation[numMassBins];    /* useful to determine Z(t - lifetime) for solving the equations */
  double numStarsAtFormation[numMassBins];
  int lower = 0, upper = 0;               /* identification of lower and upper masses for integrals */
  int index_low = 0, index_up = 0;        /* index corresponding to the specified initial metallicity */
  int *indexes = NULL;                   /* useful to avoid double computation inside nested loops */
  
  double metalYield       = 0.;   /* total metal mass released */
  double metalYieldO      = 0.;   /* total metal mass released for OXYGEN */
  double metalYieldFe     = 0.;   /* total metal mass released for IRON */
  double gasYield         = 0.;   /* total gas mass released */
  double dustYield        = 0.;   /* total dust mass released */
  
  double metalYieldRate_AGB_SNII        = 0.;   /* rate of metals produced (see GCE  equation) for AGB + SNII */
  double metalYieldRate_SNIa            = 0.;   /* rate of metals produced (see GCE  equation) for SNIa */
  double metalYieldRateO_AGB_SNII       = 0.;   /* rate of metals produced (see GCE  equation) for AGB + SNII for OXYGEN */
  double metalYieldRateO_SNIa           = 0.;   /* rate of metals produced (see GCE  equation) for SNIa       for OXYGEN */
  double metalYieldRateFe_AGB_SNII      = 0.;   /* rate of metals produced (see GCE  equation) for AGB + SNII for IRON */
  double metalYieldRateFe_SNIa          = 0.;   /* rate of metals produced (see GCE  equation) for SNIa       for IRON */

  double gasYieldRate              = 0.;   /* rate of gas    restored (see GCE  equation) for AGB + SNII */
  double gasYieldRate_SNIa         = 0.;   /* rate of gas    restored (see GCE  equation) for SNIa */
  
  double dustYieldRate             = 0.;
  double dustProductionRate_SNII   = 0.;
  double dustDestructionRate_SNII  = 0.;
//   double dustAstrationRate         = 0.;
  double dustGraingrowthRate       = 0.;
  double metallicity               = 0.;   /* metallicity in solar units */
  double dustToGasRatio            = 0.;   /* dust-to-gas ratio */
  double MSNblastwave              = 6.8e3; /* material that is accelerated to 100 km/s by SN blast wave */

  int endSnap = simParam->endSnap;
  double timeSnapsInGyr[endSnap];
  double stellarmasshistory[endSnap];
  double metalmasshistory[endSnap];

  double A    = 0.028;
  double f316 = 0.0258;
  double k    = 2.846;
  
  double *masses = simParam->metalTables_mass;
  double *metallicities = simParam->metalTables_metallicity;
  double *lifetime = simParam->metalTables_lifetime;
  double *imf = simParam->metalTables_imf;
  double *metalYieldsTable = simParam->metalTables_metalYields;
  double *metalYieldsTable_SNIa = simParam->metalTables_metalYields_SNIa;
  double *metalYieldsOTable = simParam->metalTables_metalYieldsO;
  double *metalYieldsOTable_SNIa = simParam->metalTables_metalYieldsO_SNIa;
  double *metalYieldsFeTable = simParam->metalTables_metalYieldsFe;
  double *metalYieldsFeTable_SNIa = simParam->metalTables_metalYieldsFe_SNIa;
  double *dustYieldsTable = simParam->metalTables_dustYields;
  double *gasYieldsTable = simParam->metalTables_gasYields;
  
  memset(sfrAtFormation, 0., sizeof(sfrAtFormation));
  memset(metallicityAtFormation, 0., sizeof(metallicityAtFormation));
  memset(numStarsAtFormation, 0., sizeof(numStarsAtFormation));
  memset(timeSnapsInGyr, 0., sizeof(timeSnapsInGyr));
  memset(stellarmasshistory, 0., sizeof(stellarmasshistory));
  memset(metalmasshistory, 0., sizeof(metalmasshistory));
  
  /***********************************/
  /* Compute the time and SFR arrays */
  /*  units for timeSnapsInGyr: Gyrs */
  /* conversion factor: 3.170577e-17 */
  /*   3.170577e-17 = 1/3.154e7/1e9  */
  /***********************************/
  for(int i=0; i<endSnap; i++)
  {
    timeSnapsInGyr[i] = secInGyr * simParam->timeSnaps[i];
    if(i<thisGal->snapnumber)
    {
      stellarmasshistory[i] = YrInSec * thisGal->stellarmasshistory[i] * simParam->invdeltat[i]; 
      metalmasshistory[i] = thisGal->metalmasshistory[i];
    }
  }
  // the current SFR goes to the last bin
  stellarmasshistory[thisGal->snapnumber] = sfr;
  time = timeSnapsInGyr[thisGal->snapnumber];

  /***************************************************************/
  /* Creation of array psi(t - lifetime), called sfrAtFormation       */
  /* Creation of array Z(t - lifetime), called metallicityAtFormation */
  /***************************************************************/
  for(int i=0; i<numMassBins; i++) 
  {
    if((time - lifetime[i]) > 0)
    {
      // control to avoid values outside interpolation range
      if((time - lifetime[i]) < timeSnapsInGyr[0] )
      {
        sfrAtFormation[i] = 0.;
        metallicityAtFormation[i] = 0.; 
      }
      else if((time - lifetime[i]) > timeSnapsInGyr[endSnap-1])
      {
        sfrAtFormation[i] = 0.;
        metallicityAtFormation[i] = 0.;
      }
      else
      {
        sfrAtFormation[i] = linear_interp(timeSnapsInGyr, stellarmasshistory,  endSnap, time-lifetime[i]);
        metallicityAtFormation[i] = linear_interp(timeSnapsInGyr, metalmasshistory, endSnap, time-lifetime[i]);
      }
    }
    numStarsAtFormation[i] = sfrAtFormation[i] * imf[i];
  }

  /***********************************************************************/
  /* Computation of the AGB + SNII part integral for METAL mass released */
  /*    Upper mass is 50 because after this limit yields are all zero    */
  /***********************************************************************/
  lower      = where(masses, numMassBins, 0.85);
  upper      = where(masses, numMassBins, 50.);
  // Finding the index corresponding to the specified initial metallicity
  index_low  = where(metallicities, numMetallicityBins, metallicityAtFormation[lower]);
  index_up   = where(metallicities, numMetallicityBins, metallicityAtFormation[upper]);
  indexes    = compute_indexes(lower, upper, numMetallicityBins, metallicities, metallicityAtFormation);

  metalYieldRate_AGB_SNII   = compute_yield(metalYieldsTable, numStarsAtFormation, numMassBins, indexes, index_low, index_up, lower, upper, masses[1] - masses[0]);
  metalYieldRateO_AGB_SNII  = compute_yield(metalYieldsOTable, numStarsAtFormation, numMassBins, indexes, index_low, index_up, lower, upper, masses[1] - masses[0]);
  metalYieldRateFe_AGB_SNII = compute_yield(metalYieldsFeTable, numStarsAtFormation, numMassBins, indexes, index_low, index_up, lower, upper, masses[1] - masses[0]);
  
  /*********************************************************************/
  /* Computation of the AGB + SNII part integral for GAS masses released */
  /*********************************************************************/
  gasYieldRate = compute_yield(gasYieldsTable, numStarsAtFormation, numMassBins, indexes, index_low, index_up, lower, upper, masses[1] - masses[0]);
    
  /************************************************************/
  /* Computation of all dust production and destruction rates */
  /************************************************************/
  if(thisGal->Mgas > 0.)
  {
    metallicity    = thisGal->Mmetal[0] / thisGal->Mgas;
    dustToGasRatio = thisGal->Mdust  / thisGal->Mgas; 
  }
  
  dustProductionRate_SNII  = simParam->dust_yield_SNII; /* dust production from SNII */
  dustDestructionRate_SNII = (1. - simParam->fracColdGas) * dustToGasRatio * simParam->dust_destrEfficiency * MSNblastwave;     /* dust destruction by SN blastwaves */
//   dustAstrationRate        = dustToGasRatio * sfr; /* dust astration */
  dustGraingrowthRate      = (metallicity - dustToGasRatio) * simParam->fracColdGas * thisGal->Mdust / (simParam->dust_timescaleGrainGrowth * Zsun);    /* grain growth through accretion of heavy elements in the ISM */
  dustYieldRate           += dustGraingrowthRate;
  /* dust ejection, dust astration is dealt with in produce_and_eject_metals in evol_gal.c */
  
  /********************************************/
  /* Computation of the AGB integral for DUST */
  /********************************************/
  dustYieldRate += compute_yield(dustYieldsTable, numStarsAtFormation, numMassBins, indexes, index_low, index_up, lower, upper, masses[1] - masses[0]);
    
  /***********************************************************************/
  /* Computation of DUST integral for SNII (assuming simParam->yield_dust_SNII for the yield) */
  /***********************************************************************/
  lower = where(masses, numMassBins, 8);
  upper = where(masses, numMassBins, 40);
  dustYieldRate += compute(dustProductionRate_SNII - dustDestructionRate_SNII, numStarsAtFormation, lower, upper, masses[1] - masses[0]);
  
  /*****************************************************************/
  /* Computation of the SNIa part integral for METAL mass released */
  /*****************************************************************/
  lower = where(masses, numMassBins, 0.85);
  upper = where(masses, numMassBins, 8.);
  // Computing the integral with the trapezoidal rule
  // Be careful, in this case the binning is not uniform!
  for(int i=(lower+1); i<=upper; i++)
  {
    //assuming a power law DTD (Delayed Time Distribution)
    metalYieldRate_SNIa += 0.5 * 0.15242 * ((sfrAtFormation[i] * pow(lifetime[i],-1.12)) + (sfrAtFormation[i-1] * pow(lifetime[i-1],-1.12))) * (lifetime[i] - lifetime[i-1]);
  }
  // be careful to minus sign, because lifetime is a decreasing monotonic function!
  metalYieldRateO_SNIa  = -A * f316 * k * metalYieldsOTable_SNIa[0]  * metalYieldRate_SNIa;
  metalYieldRateFe_SNIa = -A * f316 * k * metalYieldsFeTable_SNIa[0] * metalYieldRate_SNIa;
  metalYieldRate_SNIa   = -A * f316 * k * metalYieldsTable_SNIa[0]   * metalYieldRate_SNIa;
  
  /***************************************************************/
  /* Computation of the SNIa part integral for GAS mass released */
  /***************************************************************/
  // be careful because m_gas^SNIa is 0.01 solar masses
  // look at Section 2.2 of A&A 578, A87 (2015), hence:
  gasYieldRate_SNIa = (metalYieldRate_SNIa / metalYieldsTable_SNIa[0]) * 0.01;

  /**************************************************************/
  /* Total mass of METALS and GAS, multiplied by deltat because */
  /*   metalYield and gasYield are rates (see GCE equation)    */
  /*        units for deltat: years, hence 3.170577e-8          */
  /**************************************************************/
  metalYield   = secInYr * deltat * (metalYieldRate_AGB_SNII   + metalYieldRate_SNIa);
  metalYieldO  = secInYr * deltat * (metalYieldRateO_AGB_SNII  + metalYieldRateO_SNIa);
  metalYieldFe = secInYr * deltat * (metalYieldRateFe_AGB_SNII + metalYieldRateFe_SNIa);
  gasYield     = secInYr * deltat * (gasYieldRate + gasYieldRate_SNIa);
  dustYield    = secInYr * deltat * dustYieldRate;
  
  // Returned GAS mass
  *G_t = gasYield;

  // Returned METAL mass
  if (thisGal->Mgas == 0)
  {
    *e_Z    = 0.;
    *e_Z_O  = 0.;
    *e_Z_Fe = 0.; 
    *y_d    = 0.;
  }
  else 
  {
    *e_Z    = metalYield;
    *e_Z_O  = metalYieldO;
    *e_Z_Fe = metalYieldFe;
    *y_d    = dustYield;
  }
  
  free(indexes);
}
#endif
