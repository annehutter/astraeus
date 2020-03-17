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

float *getProperty(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  float invHubble = 1./simParam->hubble_h;
  float boxsize = (float)(simParam->boxsize);
  int gridsize = simParam->gridsize;
  float factor = (float) gridsize / boxsize;
  int x=0, y=0, z=0;
  
  int listToPutIn = 0;
  int offset = 0;
  outgal_t thisGal;
  float deltaTime = 1.;
  float YrInSec = 3.1536e7;
    
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);    
  float *property = allocate_array_float(numGal, "property");
  if(strcmp(propertyName,"MUV") == 0)
  {
    free(property);
    property = allocate_array_float(numGal * (currSnap + 1), "property");
  }
  
  for(int tree=0; tree<numTrees; tree++)
  {
    float **tmpProperty = malloc(sizeListEquals[tree] * sizeof(float*));
    if(tmpProperty == NULL)
    {
      fprintf(stderr, "tmpProperty could not be allocated\n");
      exit(EXIT_FAILURE);
    }
    for(int group=0; group<sizeListEquals[tree]; group++)
    {
      tmpProperty[group] = malloc((currSnap + 1) * sizeof(float));
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
        
        tmpProperty[listToPutIn][thisGal.snapnumber] = x + y*gridsize + z*gridsize*gridsize;
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
      if(strcmp(propertyName,"Mb") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * (thisGal.Mgas + thisGal.Mstar);
      }
      if(strcmp(propertyName,"Mvir") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mvir;
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
    float *UVmag = NULL;
    UVmag = calc_UV_mag(simParam, numGal, currSnap + 1, property, currSnap);
    free(property);
    property = UVmag;
  }
  
  return property;
}

float *getPropertyHistory(dconfObj_t simParam, outgtree_t **treeList, int numTrees, int ***listEquals, int **index, int *sizeListEquals, char *propertyName, int currSnap, float *times)
{
  float invHubble = 1./simParam->hubble_h;
  float boxsize = (float)(simParam->boxsize);
  int gridsize = simParam->gridsize;
  float factor = (float) gridsize / boxsize;
  int x=0, y=0, z=0;
  
  int listToPutIn = 0;
  int offset = 0;
  outgal_t thisGal;
  float deltaTime = 1.;
  float YrInSec = 3.1536e7;
    
  int numGal = get_numGal_at_snap(sizeListEquals, numTrees);
  float *property = allocate_array_float(numGal * (currSnap + 1), "property");
  
  for(int tree=0; tree<numTrees; tree++)
  {
    float **tmpProperty = malloc(sizeListEquals[tree] * sizeof(float*));
    if(tmpProperty == NULL)
    {
      fprintf(stderr, "tmpProperty could not be allocated\n");
      exit(EXIT_FAILURE);
    }
    for(int group=0; group<sizeListEquals[tree]; group++)
    {
      tmpProperty[group] = malloc((currSnap + 1) * sizeof(float));
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
      if(strcmp(propertyName,"Mb") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * (thisGal.Mgas + thisGal.Mstar);
      }
      if(strcmp(propertyName,"Mvir") == 0)
      {
        tmpProperty[listToPutIn][thisGal.snapnumber] += invHubble * thisGal.Mvir;
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
    
    float *UVmag = allocate_array_float(numGal * (currSnap + 1), "UVmag");
    for(int snap=0; snap<=currSnap; snap++)
    {
      float *UVmagSnap = calc_UV_mag(simParam, numGal, currSnap + 1, property, snap);
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

