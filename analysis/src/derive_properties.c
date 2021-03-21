#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "utils.h"
#include "dconfObj.h"
#include "outgal.h"
#include "build_index_tree_walking.h"
#include "uv.h"

int get_numGal_at_snap(int *sizeListEquals, int numTrees)
{
  int sum = 0;
  for(int tree=0; tree<numTrees; tree++)
    sum += sizeListEquals[tree];
  
  return sum;
}

double *getProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  double invHubble = 1./simParam->hubble_h;
  double barFrac = simParam->omega_b/simParam->omega_m;
  double boxsize = (double)(simParam->boxsize);
  int gridsize = simParam->gridsize;
  double factor = (double) gridsize / boxsize;
  int x=0, y=0, z=0;
  
  int listToPutIn = 0;
  int offset = 0;
  outgal_t thisGal;
  float deltaTime = 1.;
  float YrInSec = 3.1536e7;
    
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);    
  double *property = allocate_array_double(numGal, "property");
  if(strcmp(propertyName,"MUV") == 0)
  {
    free(property);
    property = allocate_array_double(numGal * (currSnap + 1), "property");
  }
  
  for(int tree=0; tree<numTrees; tree++)
  {
    double **tmpProperty = malloc(sizeListEquals[tree] * sizeof(double*));
    if(tmpProperty == NULL)
    {
      fprintf(stderr, "tmpProperty could not be allocated\n");
      exit(EXIT_FAILURE);
    }
    for(int group=0; group<sizeListEquals[tree]; group++)
    {
      tmpProperty[group] = malloc((currSnap + 1) * sizeof(double));
      if(tmpProperty[group] == NULL)
      {
        fprintf(stderr, "tmpProperty[group] could not be allocated\n");
        exit(EXIT_FAILURE);
      }
      for(int snap=0; snap<currSnap+1; snap++)
        tmpProperty[group][snap] = 0.;
    }
    
    int gal = 0;
    while(treeList[tree]->galaxies[gal].snapnumber < currSnap + 1 && gal < treeList[tree]->numGal)
    {
      thisGal = treeList[tree]->galaxies[gal];
      listToPutIn = search_index(index[tree][gal], listEquals[tree]);

      if(strcmp(propertyName, "POS") == 0)
      {
        x = (int)(factor * thisGal.pos[0]);
        y = (int)(factor * thisGal.pos[1]);
        z = (int)(factor * thisGal.pos[2]);
        
        if(x >= gridsize) x = gridsize-1;
        if(y >= gridsize) y = gridsize-1;
        if(z >= gridsize) z = gridsize-1;
        
        if(x >= gridsize || y >= gridsize || z >= gridsize)
          printf("PROBLEM! x = %d  y = %d  z = %d\t rank = %d snap = %d listToPutIn = %d\t tree = %d gal = %d offset = %d\t POS = %d\n", x, y, z, simParam->thisRank, thisGal.snapnumber, listToPutIn, tree, gal, offset, x + y*gridsize + z*gridsize*gridsize);
     
        tmpProperty[listToPutIn][thisGal.snapnumber] = x + y*gridsize + z*gridsize*gridsize;
      }
      if(strcmp(propertyName, "POSX") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.pos[0];
      }
      if(strcmp(propertyName, "POSY") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.pos[1];
      }
      if(strcmp(propertyName, "POSZ") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.pos[2];
      }
      if(strcmp(propertyName, "VELX") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.vel[0];
      }
      if(strcmp(propertyName, "VELY") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.vel[1];
      }
      if(strcmp(propertyName, "VELZ") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.vel[2];
      }
      if(strcmp(propertyName, "SFR") == 0)
      {
        if(thisGal.snapnumber > 0) 
          deltaTime = times[thisGal.snapnumber] - times[thisGal.snapnumber - 1];
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.MgasIni * thisGal.feff / deltaTime * YrInSec;
      }
      if(strcmp(propertyName, "MUV") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.MgasIni * thisGal.feff;
      }
      if(strcmp(propertyName, "newMstar") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.MgasIni * thisGal.feff;
      }
      if(strcmp(propertyName,"Mstar") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mstar;
      }
      if(strcmp(propertyName,"Mgas") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mgas;
      }
      if(strcmp(propertyName,"MgasIni") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.MgasIni;
      }
      if(strcmp(propertyName,"MgasAcc") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * barFrac * (thisGal.Mvir - thisGal.Mvir_prog);
      }
      if(strcmp(propertyName,"MgasAccRate") == 0)
      {
        if(thisGal.snapnumber > 0) 
          deltaTime = times[thisGal.snapnumber] - times[thisGal.snapnumber - 1];
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * barFrac * (thisGal.Mvir - thisGal.Mvir_prog) / deltaTime * YrInSec;
      }
      if(strcmp(propertyName,"MgasMer") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * (thisGal.MgasIni - barFrac * (thisGal.Mvir - thisGal.Mvir_prog));
      }
      if(strcmp(propertyName,"Mb") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * (thisGal.Mgas + thisGal.Mstar);
      }
      if(strcmp(propertyName,"Mvir") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mvir;
      }
      if(strcmp(propertyName,"MvirProg") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mvir_prog;
      }
      if(strcmp(propertyName,"Rvir") == 0)
      {
        if(invHubble * thisGal.Rvir > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Rvir;
      }
      if(strcmp(propertyName,"Spin") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.spin;
      }
      if(strcmp(propertyName,"feff") == 0)
        if(thisGal.feff > tmpProperty[listToPutIn][thisGal.snapnumber])
        {
          tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.feff;
        }
      if(strcmp(propertyName,"fg") == 0)
      {
        if(thisGal.fg > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.fg;
      }
      if(strcmp(propertyName,"zreion") == 0)
      {
        if(1./thisGal.scalefactor-1. > thisGal.zreion && thisGal.zreion > 0.)
          printf("SOMETHING WRONG: z = %e\t zreion = %e\n", 1./thisGal.scalefactor-1., thisGal.zreion);
        if(thisGal.zreion > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.zreion;
      }
      if(strcmp(propertyName,"photHI_bg") == 0)
      {
        if(thisGal.photHI_bg > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.photHI_bg;
      }
      
      if(gal < treeList[tree]->numGal - 1)
        gal++;
      else
        break;
    }
    
    if(strcmp(propertyName,"MUV") == 0)
    {
      for(int group=offset; group<offset + sizeListEquals[tree]; group++)
      {
        for(int snap=0; snap<=currSnap; snap++)
        {
          property[group*(currSnap + 1) + snap] = tmpProperty[group-offset][snap];
        }
      }
    }
    else
    {
      for(int group=offset; group<offset + sizeListEquals[tree]; group++)
      {
        property[group] = tmpProperty[group-offset][currSnap];
      }
    }
    
    offset += sizeListEquals[tree];
    
    for(int group=0; group<sizeListEquals[tree]; group++)
      free(tmpProperty[group]);
    free(tmpProperty);
  }
  
  if(strcmp(propertyName, "MUV") == 0)
  {
    double *UVmag = NULL;
    UVmag = calc_UV_mag(simParam, numGal, currSnap + 1, property, currSnap);
    free(property);
    property = UVmag;
  }
  
  return property;
}

double *getPropertyHistory(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  double invHubble = 1./simParam->hubble_h;
  double barFrac = simParam->omega_b/simParam->omega_m;
  double boxsize = (double)(simParam->boxsize);
  int gridsize = simParam->gridsize;
  double factor = (double) gridsize / boxsize;
  int x=0, y=0, z=0;
  
  int listToPutIn = 0;
  int offset = 0;
  outgal_t thisGal;
  float deltaTime = 1.;
  float YrInSec = 3.1536e7;
    
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
  double *property = allocate_array_double(numGal * (currSnap + 1), "property");
  
  for(int tree=0; tree<numTrees; tree++)
  {
    double **tmpProperty = malloc(sizeListEquals[tree] * sizeof(double*));
    if(tmpProperty == NULL)
    {
      fprintf(stderr, "tmpProperty could not be allocated\n");
      exit(EXIT_FAILURE);
    }
    for(int group=0; group<sizeListEquals[tree]; group++)
    {
      tmpProperty[group] = malloc((currSnap + 1) * sizeof(double));
      if(tmpProperty[group] == NULL)
      {
        fprintf(stderr, "tmpProperty[group] could not be allocated\n");
        exit(EXIT_FAILURE);
      }
      for(int snap=0; snap<currSnap+1; snap++)
        tmpProperty[group][snap] = -1.;
    }
    
    int gal = 0;
    while(treeList[tree]->galaxies[gal].snapnumber < currSnap + 1 && gal < treeList[tree]->numGal)
    {
      thisGal = treeList[tree]->galaxies[gal];
      listToPutIn = search_index(index[tree][gal], listEquals[tree]);
      
      if(tmpProperty[listToPutIn][thisGal.snapnumber] == -1.)
        tmpProperty[listToPutIn][thisGal.snapnumber] = 0.;

      if(strcmp(propertyName, "POS") == 0)
      {
        x = (int)(factor * thisGal.pos[0]);
        y = (int)(factor * thisGal.pos[1]);
        z = (int)(factor * thisGal.pos[2]);

        if(x >= gridsize) x = gridsize-1;
        if(y >= gridsize) y = gridsize-1;
        if(z >= gridsize) z = gridsize-1;
        
        tmpProperty[listToPutIn][thisGal.snapnumber] = x + y*gridsize + z*gridsize*gridsize;
      }
      if(strcmp(propertyName, "POSX") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.pos[0];
      }
      if(strcmp(propertyName, "POSY") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.pos[1];
      }
      if(strcmp(propertyName, "POSZ") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.pos[2];
      }
      if(strcmp(propertyName, "VELX") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.vel[0];
      }
      if(strcmp(propertyName, "VELY") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.vel[1];
      }
      if(strcmp(propertyName, "VELZ") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.vel[2];
      }
      if(strcmp(propertyName, "SFR") == 0)
      {
        if(thisGal.snapnumber > 0) 
          deltaTime = times[thisGal.snapnumber] - times[thisGal.snapnumber - 1];
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.MgasIni * thisGal.feff / deltaTime * YrInSec;
      }
      if(strcmp(propertyName, "MUV") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.MgasIni * thisGal.feff;
      }
      if(strcmp(propertyName, "newMstar") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.MgasIni * thisGal.feff;
      }
      if(strcmp(propertyName,"Mstar") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mstar;
      }
      if(strcmp(propertyName,"Mgas") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mgas;
      }
      if(strcmp(propertyName,"MgasIni") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.MgasIni;
      }
      if(strcmp(propertyName,"MgasAcc") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * barFrac * (thisGal.Mvir - thisGal.Mvir_prog);
      }
      if(strcmp(propertyName,"MgasAccRate") == 0)
      {
        if(thisGal.snapnumber > 0) 
          deltaTime = times[thisGal.snapnumber] - times[thisGal.snapnumber - 1];
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * barFrac * (thisGal.Mvir - thisGal.Mvir_prog) / deltaTime * YrInSec;
      }
      if(strcmp(propertyName,"MgasMer") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * (thisGal.MgasIni - barFrac * (thisGal.Mvir - thisGal.Mvir_prog));
      }
      if(strcmp(propertyName,"Mb") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * (thisGal.Mgas + thisGal.Mstar);
      }
      if(strcmp(propertyName,"Rvir") == 0)
      {
        if(invHubble * thisGal.Rvir > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = invHubble * thisGal.Rvir;
      }
      if(strcmp(propertyName,"Mvir") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mvir;
      }
      if(strcmp(propertyName,"MvirProg") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mvir_prog;
      }
      if(strcmp(propertyName,"Spin") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.spin;
      }
      if(strcmp(propertyName,"feff") == 0)
        if(thisGal.feff > tmpProperty[listToPutIn][thisGal.snapnumber])
        {
          tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.feff;
        }
      if(strcmp(propertyName,"fg") == 0)
      {
        if(thisGal.fg > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.fg;
      }
      if(strcmp(propertyName,"zreion") == 0)
      {
        if(thisGal.zreion > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.zreion;
      }
      if(strcmp(propertyName,"photHI_bg") == 0)
      {
        if(thisGal.photHI_bg > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.photHI_bg;
      }
      
      if(gal < treeList[tree]->numGal - 1)
        gal++;
      else
        break;
    }
    
    for(int group=offset; group<offset + sizeListEquals[tree]; group++)
    {
      for(int snap=0; snap<=currSnap; snap++)
      {
        property[group*(currSnap + 1) + snap] = tmpProperty[group-offset][snap];
      }
    }
    offset += sizeListEquals[tree];
    
    for(int group=0; group<sizeListEquals[tree]; group++)
      free(tmpProperty[group]);
    free(tmpProperty);
  }
  
  if(strcmp(propertyName, "MUV") == 0)
  {
    for(int gal=0; gal<numGal; gal++)
    {
      for(int snap=0; snap<=currSnap; snap++)
      {
        if(property[gal*(currSnap + 1) + snap] == -1.)
          property[gal*(currSnap + 1) + snap] = 0.;
      }
    }
    
    double *UVmag = allocate_array_double(numGal * (currSnap + 1), "UVmag");
    for(int snap=0; snap<=currSnap; snap++)
    {
      double *UVmagSnap = calc_UV_mag(simParam, numGal, currSnap + 1, property, snap);
      for(int gal=0; gal<numGal; gal++)
      {
        UVmag[gal*(currSnap + 1) + snap] = UVmagSnap[gal];
      }
      free(UVmagSnap);
    }
    
    free(property);
    property = UVmag;
  }
  
  return property;
}

double *getProperty_endSnap(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  double invHubble = 1./simParam->hubble_h;
  double barFrac = simParam->omega_b/simParam->omega_m;
  double boxsize = (double)(simParam->boxsize);
  int gridsize = simParam->gridsize;
  double factor = (double) gridsize / boxsize;
  int x=0, y=0, z=0;
  
  int listToPutIn = 0;
  int offset = 0;
  outgal_t thisGal;
  outgal_t endGal;
  float deltaTime = 1.;
  float YrInSec = 3.1536e7;
    
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);    
  double *property = allocate_array_double(numGal, "property");
  if(strcmp(propertyName,"MUV") == 0)
  {
    free(property);
    property = allocate_array_double(numGal * (currSnap + 1), "property");
  }
  
  for(int tree=0; tree<numTrees; tree++)
  {
    double **tmpProperty = malloc(sizeListEquals[tree] * sizeof(double*));
    if(tmpProperty == NULL)
    {
      fprintf(stderr, "tmpProperty could not be allocated\n");
      exit(EXIT_FAILURE);
    }
    for(int group=0; group<sizeListEquals[tree]; group++)
    {
      tmpProperty[group] = malloc((currSnap + 1) * sizeof(double));
      if(tmpProperty[group] == NULL)
      {
        fprintf(stderr, "tmpProperty[group] could not be allocated\n");
        exit(EXIT_FAILURE);
      }
      for(int snap=0; snap<currSnap+1; snap++)
        tmpProperty[group][snap] = 0.;
    }
    
    endGal = treeList[tree]->galaxies[treeList[tree]->numGal-1];
    
    int gal = 0;
    while(treeList[tree]->galaxies[gal].snapnumber < currSnap + 1 && gal < treeList[tree]->numGal)
    {
      thisGal = treeList[tree]->galaxies[gal];
      listToPutIn = search_index(index[tree][gal], listEquals[tree]);

      if(strcmp(propertyName, "POS") == 0)
      {
        x = (int)(factor * endGal.pos[0]);
        y = (int)(factor * endGal.pos[1]);
        z = (int)(factor * endGal.pos[2]);
        
        if(x >= gridsize) x = gridsize-1;
        if(y >= gridsize) y = gridsize-1;
        if(z >= gridsize) z = gridsize-1;
        
        if(x >= gridsize || y >= gridsize || z >= gridsize)
          printf("PROBLEM! x = %d  y = %d  z = %d\t rank = %d snap = %d listToPutIn = %d\t tree = %d gal = %d offset = %d\t POS = %d\n", x, y, z, simParam->thisRank, endGal.snapnumber, listToPutIn, tree, gal, offset, x + y*gridsize + z*gridsize*gridsize);
     
        tmpProperty[listToPutIn][thisGal.snapnumber] = x + y*gridsize + z*gridsize*gridsize;
      }
      if(strcmp(propertyName, "POSX") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.pos[0];
      }
      if(strcmp(propertyName, "POSY") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.pos[1];
      }
      if(strcmp(propertyName, "POSZ") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.pos[2];
      }
      if(strcmp(propertyName, "VELX") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.vel[0];
      }
      if(strcmp(propertyName, "VELY") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.vel[1];
      }
      if(strcmp(propertyName, "VELZ") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.vel[2];
      }
      if(strcmp(propertyName, "SFR") == 0)
      {
        if(thisGal.snapnumber > 0) 
          deltaTime = times[endGal.snapnumber] - times[endGal.snapnumber - 1];
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * endGal.MgasIni * endGal.feff / deltaTime * YrInSec;
      }
      if(strcmp(propertyName, "newMstar") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * endGal.MgasIni * endGal.feff;
      }
      if(strcmp(propertyName,"Mstar") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * endGal.Mstar;
      }
      if(strcmp(propertyName,"Mgas") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * endGal.Mgas;
      }
      if(strcmp(propertyName,"MgasIni") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * endGal.MgasIni;
      }
      if(strcmp(propertyName,"MgasAcc") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * barFrac * (endGal.Mvir - endGal.Mvir_prog);
      }
      if(strcmp(propertyName,"MgasAccRate") == 0)
      {
        if(thisGal.snapnumber > 0) 
          deltaTime = times[endGal.snapnumber] - times[endGal.snapnumber - 1];
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * barFrac * (endGal.Mvir - endGal.Mvir_prog) / deltaTime * YrInSec;
      }
      if(strcmp(propertyName,"MgasMer") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * (endGal.MgasIni - barFrac * (endGal.Mvir - endGal.Mvir_prog));
      }
      if(strcmp(propertyName,"Mb") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * (endGal.Mgas + endGal.Mstar);
      }
      if(strcmp(propertyName,"Mvir") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * endGal.Mvir;
      }
      if(strcmp(propertyName,"MvirProg") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * endGal.Mvir_prog;
      }
      if(strcmp(propertyName,"Rvir") == 0)
      {
        if(invHubble * endGal.Rvir > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = invHubble * endGal.Rvir;
      }
      if(strcmp(propertyName,"Spin") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] = thisGal.spin;
      }
      if(strcmp(propertyName,"feff") == 0)
        if(endGal.feff > tmpProperty[listToPutIn][thisGal.snapnumber])
        {
          tmpProperty[listToPutIn][thisGal.snapnumber] = endGal.feff;
        }
      if(strcmp(propertyName,"fg") == 0)
      {
        if(endGal.fg > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = endGal.fg;
      }
      if(strcmp(propertyName,"zreion") == 0)
      {
        if(1./endGal.scalefactor-1. > endGal.zreion && endGal.zreion > 0.)
          printf("SOMETHING WRONG: z = %e\t zreion = %e\n", 1./endGal.scalefactor-1., endGal.zreion);
        if(endGal.zreion > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = endGal.zreion;
      }
      if(strcmp(propertyName,"photHI_bg") == 0)
      {
        if(endGal.photHI_bg > tmpProperty[listToPutIn][thisGal.snapnumber])
          tmpProperty[listToPutIn][thisGal.snapnumber] = endGal.photHI_bg;
      }
      
      if(gal < treeList[tree]->numGal - 1)
        gal++;
      else
        break;
    }
    
    for(int group=offset; group<offset + sizeListEquals[tree]; group++)
    {
      property[group] = tmpProperty[group-offset][currSnap];
    }
    
    offset += sizeListEquals[tree];
    
    for(int group=0; group<sizeListEquals[tree]; group++)
      free(tmpProperty[group]);
    free(tmpProperty);
  }
  
  return property;
}
