#ifndef OUTGAL_H
#define OUTGAL_H

typedef struct{
  float scalefactor;
  
  float pos[3];
  float vel[3];
  float Mvir;
  float Mvir_prog;
  float Rvir;
  float velDisp;
  float velMax;
  float spin;
  float scalefactorLastMajorMerger;
  float MgasIni;
  float Mgas;
  float Mstar;
  
  float feff;
  float fg;
  float zreion;
  float photHI_bg;
  float stellarmasshistory[OUTPUTSNAPNUMBER];
}outgalsnap_t;


typedef struct{
  int32_t localID;
  int32_t localDescID;
  int32_t numProg;
  int32_t snapnumber;
  
  float scalefactor;
  float pos[3];
  float vel[3];
  float Mvir;
  float Mvir_prog;
  float Rvir;
  float velDisp;
  float velMax;
  float spin;
  float scalefactorLastMajorMerger;
  float MgasIni;
  float Mgas;
  float Mstar;
  
  float feff;
  float fg;
  float zreion;
  float photHI_bg;
}outgal_t; 

typedef struct{
  int32_t numGal;
  outgal_t *outgalaxies;
}outgtree_t;


/*-------------------------------------------------------*/
/* FUNCTIONS ON OUTGALSNAP_T                             */
/*-------------------------------------------------------*/

outgalsnap_t *initOutGalSnap();
void initOutGalSnapWithoutAllocation(outgalsnap_t *newOutGal);
outgalsnap_t *initOutGalList(int numGal);
void deallocate_outgalsnap(outgalsnap_t *thisOutGal);

/*-------------------------------------------------------*/
/* FUNCTIONS ON OUTGAL_T                                 */
/*-------------------------------------------------------*/

outgal_t *initOutGalTree();
void deallocate_outgaltree(outgal_t *thisOutGal);

/*-------------------------------------------------------*/
/* FUNCTIONS ON GTREE_T                                  */
/*-------------------------------------------------------*/

outgtree_t *initOutGtree(int numGal);
void deallocate_outgtree(outgtree_t *thisOutGtree);

#endif
