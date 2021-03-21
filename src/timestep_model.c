#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>
#include <math.h>

#include "dconfObj.h"
#include "timestep_model.h"

float get_corrFactorTimeStep(dconfObj_t simParam, int snap)
{
  float corrFactorTimeStep = 1.;
  
  if(simParam->dt_model == 1)
  {
    corrFactorTimeStep = simParam->dt_rescaleFactor;
  }
  if(simParam->dt_model > 1 && snap > 0)
  {
    corrFactorTimeStep = (simParam->timeSnaps[snap] - simParam->timeSnaps[snap-1]) * 3.170979e-14 / simParam->dt_deltaTimeInMyr;
  }
  
  return corrFactorTimeStep;
}