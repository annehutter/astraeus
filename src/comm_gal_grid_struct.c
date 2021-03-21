#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <complex.h>

#include "utils.h"
#include "comm_gal_grid_struct.h"

commGalGrid_t *initCommGalGrid(int32_t numGal)
{
  commGalGrid_t *thisCommGalGrid = malloc(sizeof(commGalGrid_t));
  
  if(thisCommGalGrid == NULL)
  {
    fprintf(stderr, "Could not allocate commGalGrid_t\n");
    exit(EXIT_FAILURE);
  }
  
  thisCommGalGrid->numGal = numGal;
  thisCommGalGrid->numGalWritten = 0;
  thisCommGalGrid->zreion = allocate_array_float(numGal, "zreion");
  thisCommGalGrid->photHI = allocate_array_float(numGal, "photHI");
  thisCommGalGrid->pos = allocate_array_int32_t(numGal, "pos");

  return thisCommGalGrid;
}

commGalGrid_t *initCommGalGrid_pos(int32_t numGal)
{
  commGalGrid_t *thisCommGalGrid = malloc(sizeof(commGalGrid_t));
  
  if(thisCommGalGrid == NULL)
  {
    fprintf(stderr, "Could not allocate commGalGrid_t\n");
    exit(EXIT_FAILURE);
  }
  
  thisCommGalGrid->numGal = numGal;
  thisCommGalGrid->numGalWritten = 0;
  thisCommGalGrid->zreion =  NULL;
  thisCommGalGrid->photHI = NULL;
  thisCommGalGrid->pos = allocate_array_int32_t(numGal, "pos");

  return thisCommGalGrid;
}

commGalGrid_t *initCommGalGrid_zreion_photHI(int32_t numGal)
{
  commGalGrid_t *thisCommGalGrid = malloc(sizeof(commGalGrid_t));
  
  if(thisCommGalGrid == NULL)
  {
    fprintf(stderr, "Could not allocate commGalGrid_t\n");
    exit(EXIT_FAILURE);
  }
  
  thisCommGalGrid->numGal = numGal;
  thisCommGalGrid->numGalWritten = 0;
  thisCommGalGrid->zreion = allocate_array_float(numGal, "zreion");
  thisCommGalGrid->photHI = allocate_array_float(numGal, "photHI");
  thisCommGalGrid->pos = NULL;

  return thisCommGalGrid;
}

void allocCommGalGrid_zreion_photHI(commGalGrid_t **thisCommGalGrid, int32_t numGal)
{
  (*thisCommGalGrid)->numGal = numGal;
  (*thisCommGalGrid)->zreion = allocate_array_float(numGal, "zreion");
  (*thisCommGalGrid)->photHI = allocate_array_float(numGal, "photHI");
}

void reallocCommGalGrid(commGalGrid_t **thisCommGalGrid, int32_t numGal)
{  
  (*thisCommGalGrid)->numGal = numGal;
  (*thisCommGalGrid)->zreion = realloc((*thisCommGalGrid)->zreion, sizeof(float) * numGal);
  (*thisCommGalGrid)->photHI = realloc((*thisCommGalGrid)->photHI, sizeof(float) * numGal);
  (*thisCommGalGrid)->pos = realloc((*thisCommGalGrid)->pos, sizeof(int32_t) * numGal);
}

void reallocCommGalGrid_pos(commGalGrid_t **thisCommGalGrid, int32_t numGal)
{  
  (*thisCommGalGrid)->numGal = numGal;
  (*thisCommGalGrid)->pos = realloc((*thisCommGalGrid)->pos, sizeof(int32_t) * numGal);
}

void reallocCommGalGrid_zreion_photHI(commGalGrid_t **thisCommGalGrid, int32_t numGal)
{  
  (*thisCommGalGrid)->numGal = numGal;
  (*thisCommGalGrid)->zreion = realloc((*thisCommGalGrid)->zreion, sizeof(float) * numGal);
  (*thisCommGalGrid)->photHI = realloc((*thisCommGalGrid)->photHI, sizeof(float) * numGal);
}


void deallocate_commGalGrid(commGalGrid_t *thisCommGalGrid)
{
  if(thisCommGalGrid->zreion != NULL) free(thisCommGalGrid->zreion);
  if(thisCommGalGrid->photHI != NULL) free(thisCommGalGrid->photHI);
  if(thisCommGalGrid->pos != NULL) free(thisCommGalGrid->pos);

  free(thisCommGalGrid);
}


