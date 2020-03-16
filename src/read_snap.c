#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "utils.h"
#include "gal_gtree.h"
#include "read_snap.h"

outgalsnap_t *read_snap(int snap, int nb_gal, char *baseName)
{
    outgalsnap_t *galList = NULL;
    outgalsnap_t thisGal;
    FILE *f;
    int count = 0;	
    
    char infileEnding[10];
    char *thisInfileName = NULL;
    
    galList = malloc(sizeof(outgalsnap_t) * nb_gal);
    sprintf(infileEnding, "tmp_Galsnap_%d.dat", snap);
    thisInfileName = concat(baseName, infileEnding);
    f = fopen(thisInfileName, "r"); 
    if(f == NULL)
        printf("Error opening file : %s\n", thisInfileName);
      
    while(fread(&thisGal, sizeof(outgalsnap_t), 1, f))
    {
        galList[count] = thisGal;
        count++;
    }
    free(thisInfileName);
    fclose(f);
    return galList;
}

outgalsnap_t **read_snaps(int minSnap, int maxSnap, int *nb_gals)
{
    outgalsnap_t **snapList = NULL;
    snapList = malloc((maxSnap-minSnap+1) * sizeof(outgalsnap_t*));
    for(int snap=minSnap; snap<=maxSnap; snap++)
    {
        snapList[snap-minSnap] = read_snap(snap, nb_gals[snap], "/data/users/legrand/test/");
    }
    return snapList;
}

int *get_nblines()
{
    FILE *f_nbline = fopen("/data/users/legrand/test/nbline.dat", "r");
    static int nb_line[100] = {0}; 
    fread(nb_line, sizeof(int), 100, f_nbline);
    
    fclose(f_nbline);
    return nb_line;
}
