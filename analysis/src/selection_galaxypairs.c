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
#include "read_trees.h"
#include "zas_lists.h"

#include "build_index_tree_walking.h"
#include "derive_properties.h"
#include "statistics.h"
#include "operations_on_properties.h"
#include "statistics_galaxypairs.h"

#include "selection_galaxypairs.h"

void select_galaxypairs(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *selectionPropertyName, char *selectionProperty2Name, char *propertyWithHistoryName, double minSelectionProperty, double maxSelectionProperty, double minSelectionProperty2, double maxSelectionProperty2, double maxDistance, char *filename, int thisRank, int size)
{
  int numGal = 0;      
  double *posx = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "POSX", currSnap, simParam->times, &numGal);
  double *posy = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "POSY", currSnap, simParam->times, &numGal);
  double *posz = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "POSZ", currSnap, simParam->times, &numGal);
  double *Mvir = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Mvir", currSnap, simParam->times, &numGal);
  double *MgasIni = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "MgasIni", currSnap, simParam->times, &numGal);
  double *Mgas = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Mgas", currSnap, simParam->times, &numGal);
  double *Mstar = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Mstar", currSnap, simParam->times, &numGal);
  double *density = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "DENS", currSnap, simParam->times, &numGal);
  double *XHII = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "XHII", currSnap, simParam->times, &numGal);

  double *selectionProperty = NULL;
  if(strcmp(selectionPropertyName, "Mvir") == 0) selectionProperty = Mvir;
  else if(strcmp(selectionPropertyName, "MgasIni") == 0) selectionProperty = MgasIni;
  else if(strcmp(selectionPropertyName, "Mgas") == 0) selectionProperty = Mgas;
  else if(strcmp(selectionPropertyName, "Mstar") == 0) selectionProperty = Mstar;
  else if(strcmp(selectionPropertyName, "density") == 0) selectionProperty = density;
  else if(strcmp(selectionPropertyName, "XHII") == 0) selectionProperty = XHII;
  else selectionProperty = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, selectionPropertyName, currSnap, simParam->times, &numGal);
  
  double *selectionProperty2 = NULL;
  if(strcmp(selectionProperty2Name, "Mvir") == 0) selectionProperty2 = Mvir;
  else if(strcmp(selectionProperty2Name, "MgasIni") == 0) selectionProperty2 = MgasIni;
  else if(strcmp(selectionProperty2Name, "Mgas") == 0) selectionProperty2 = Mgas;
  else if(strcmp(selectionProperty2Name, "Mstar") == 0) selectionProperty2 = Mstar;
  else if(strcmp(selectionProperty2Name, "density") == 0) selectionProperty2 = density;
  else if(strcmp(selectionProperty2Name, "XHII") == 0) selectionProperty2 = XHII;
  else selectionProperty2 = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, selectionProperty2Name, currSnap, simParam->times, &numGal);
  
  double *propertyWithHistory = getThisPropertyHistory(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, propertyWithHistoryName, currSnap, simParam->times, &numGal);
  
  statistics_on_galaxypairs(numGal, selectionProperty, minSelectionProperty, maxSelectionProperty, selectionProperty2, minSelectionProperty2, maxSelectionProperty2, maxDistance, currSnap+1, simParam->redshifts, currSnap, posx, posy, posz, Mvir, MgasIni, Mgas, Mstar, density, XHII, propertyWithHistory, thisRank, size, filename);

  free(posx);
  free(posy);
  free(posz);
  
  free(Mvir);
  free(MgasIni);
  free(Mgas);
  free(Mstar);
  free(density);
  free(XHII);
  
  if(strcmp(selectionPropertyName, "Mvir") != 0 && strcmp(selectionPropertyName, "MgasIni") != 0 &&           strcmp(selectionPropertyName, "Mgas") != 0 && strcmp(selectionPropertyName, "Mstar") != 0 && strcmp(selectionPropertyName, "density") != 0 && strcmp(selectionPropertyName, "XHII") != 0) free(selectionProperty);
  if(strcmp(selectionProperty2Name, "Mvir") != 0 && strcmp(selectionProperty2Name, "MgasIni") != 0 &&           strcmp(selectionProperty2Name, "Mgas") != 0 && strcmp(selectionProperty2Name, "Mstar") != 0 && strcmp(selectionProperty2Name, "density") != 0 && strcmp(selectionProperty2Name, "XHII") != 0) free(selectionProperty2);
  free(propertyWithHistory);
}

char *create_filename_galaxypairs(char *directory, char *selectionPropertyName, char *selectionProperty2Name, double minSelectionProperty, double maxSelectionProperty, double minSelectionProperty2, double maxSelectionProperty2, float redshift)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);
  
  char selectionPropertyRangeName[24];
  sprintf(selectionPropertyRangeName, "_%3.2e-%3.2e", minSelectionProperty, maxSelectionProperty);
  char selectionProperty2RangeName[24];
  sprintf(selectionProperty2RangeName, "_%3.2e-%3.2e", minSelectionProperty2, maxSelectionProperty2);
  
  char *filename = NULL;
  filename = concat_strings(8, directory, "/galaxypairs_", selectionPropertyName, selectionPropertyRangeName, "_", selectionProperty2Name, selectionProperty2RangeName, redshift_name);
  
  return filename;
}