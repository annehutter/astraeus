#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <sys/stat.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "dconfObj.h"
#include "outgal.h"
#include "operations_on_properties.h"
#include "build_index_tree_walking.h"

#include "output_lists.h"

void write_galaxies_to_textfile_all_snaps(dconfObj_t simParam, outgtree_t **theseTrees, int numTrees, int **index, int ***listEquals, int *sizeListEquals, int **listMerged, int thisRank)
{    
  int numSnaps = simParam->numSnaps;

  int currSnap = 0, prevSnap = -1;
  for(int snap=0; snap<numSnaps; snap++)
  { 
    if(thisRank == 0) printf("Writing galaxies to text file: snap = %d\n", snap);
    
    currSnap = snap;
    get_listEquals(theseTrees, numTrees, index, listEquals, sizeListEquals, listMerged, prevSnap, currSnap);
    prevSnap = currSnap;
    
    char *filename = create_filename_output_lists(simParam->outputDir, simParam->redshifts[currSnap], thisRank);
    write_galaxies_to_textfile(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, currSnap, filename, thisRank);
    free(filename);
  }
}

void write_galaxies_to_textfile(dconfObj_t simParam, outgtree_t **theseTrees, int32_t numTrees, int ***listEquals, int **index, int *sizeListEquals, int currSnap, char *filename, int thisRank)
{
  int numGal = 0;
  
  double *Mvir = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Mvir", currSnap, simParam->times, &numGal);
  double *dens = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "DENS", currSnap, simParam->times, &numGal);
  double *posx = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "POSX", currSnap, simParam->times, &numGal);
  double *posy = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "POSY", currSnap, simParam->times, &numGal);
  double *posz = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "POSZ", currSnap, simParam->times, &numGal);
  double *velx = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "VELX", currSnap, simParam->times, &numGal);
  double *vely = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "VELY", currSnap, simParam->times, &numGal);
  double *velz = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "VELZ", currSnap, simParam->times, &numGal);
  double *spin = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Spin", currSnap, simParam->times, &numGal);
  double *Mstar = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Mstar", currSnap, simParam->times, &numGal);
  double *newMstar = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "newMstar", currSnap, simParam->times, &numGal);
  double *Mgas = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Mgas", currSnap, simParam->times, &numGal);
  double *MgasIni = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "MgasIni", currSnap, simParam->times, &numGal);
  double *SFR = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "SFR", currSnap, simParam->times, &numGal);
  double *zreion = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "zreion", currSnap, simParam->times, &numGal);
  double *XHI = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "XHI", currSnap, simParam->times, &numGal);
  double *MUV = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "MUV", currSnap, simParam->times, &numGal);
  double *Nion = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "Nion", currSnap, simParam->times, &numGal);
  double *fesc = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "fesc", currSnap, simParam->times, &numGal);
  double *fescFej = getThisProperty(simParam, theseTrees, numTrees, listEquals, index, sizeListEquals, "fescFej", currSnap, simParam->times, &numGal);

  FILE *f = NULL;
  f = fopen(filename, "w");
  if(f == NULL)
  {
    fprintf(stderr, "Could not open output file %s\n", filename);
    exit(EXIT_FAILURE);
  }
  
  if(thisRank == 0) fprintf(f, "Mvir [Msun]\t rho/<rho>\t x [h^-1 Mpc]\t y [h^-1 Mpc]\t z [h^-1 Mpc]\t vx [km/s phys]\t vy [km/s phys]\t vz [km/s phys]\t Spin\t Mstar [Msun]\t newMstar [Msun]\t Mgas [Msun]\t MgasIni [Msun]\t SFR [Msun/yr]\t zreion \t XHI\t MUV\t Nion [1.e50 s^-1]\t fesc\t fescFej\n");
  
  for(int gal=0; gal<numGal-1; gal++)
  {
    fprintf(f, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", Mvir[gal], dens[gal], posx[gal], posy[gal], posz[gal], velx[gal], vely[gal], velz[gal], spin[gal], Mstar[gal], newMstar[gal], Mgas[gal], MgasIni[gal], SFR[gal], zreion[gal], XHI[gal], MUV[gal], Nion[gal], fesc[gal], fescFej[gal]);
  }
  fprintf(f, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", Mvir[numGal-1], dens[numGal-1],  posx[numGal-1], posy[numGal-1], posz[numGal-1], velx[numGal-1], vely[numGal-1], velz[numGal-1], spin[numGal-1], Mstar[numGal-1], newMstar[numGal-1], Mgas[numGal-1], MgasIni[numGal-1], SFR[numGal-1], zreion[numGal-1], XHI[numGal-1], MUV[numGal-1], Nion[numGal-1], fesc[numGal-1], fescFej[numGal-1]);

  fclose(f);
  
  free(Mvir);
  free(dens);
  free(posx);
  free(posy);
  free(posz);
  free(velx);
  free(vely);
  free(velz);
  free(spin);
  free(Mstar);
  free(newMstar);
  free(Mgas);
  free(MgasIni);
  free(SFR);
  free(zreion);
  free(XHI);
  free(MUV);
  free(Nion);
  free(fesc);
  free(fescFej);
}

char *create_filename_output_lists(char *directory, float redshift, int thisRank)
{
  char redshift_name[12];
  sprintf(redshift_name, "_z%3.2f", redshift);
  
  char rank_name[12];
  sprintf(rank_name, "_%d", thisRank);
  
  char *subdirectory = concat_strings(2, directory, "/output_lists");
  if (mkdir(subdirectory, 0777) == -1) 
    fprintf(stderr, "Could not create directory %s\n", subdirectory);

  char *filename = concat_strings(5, subdirectory, "/galaxies_propertyLists", redshift_name, rank_name, ".txt");
  
  free(subdirectory);
  
  return filename;
}