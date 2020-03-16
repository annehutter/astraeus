#ifndef GAL_GTREE_H
#define GAL_GTREE_H

typedef struct{
  int32_t localID;
  int32_t localDescID;
  int32_t numProg;

  int32_t snapnumber;
  float scalefactor;
  float descScalefactor;
  
  float pos[3];
  float vel[3];
  float Mvir;
  float Mvir_prog;
  float Rvir;
  float velDisp;
  float velMax;
  float spin;
  float scalefactorLastMajorMerger;
  float Mgas;
  float MgasIni;
  float Mstar;
  
  float feff;
  float fg;
  float photHI_bg;
  float zreion;
  
  float *stellarmasshistory;
}gal_t;

typedef struct{
  int32_t numGal;
  int32_t walker;
  gal_t *galaxies;
}gtree_t;

/*-------------------------------------------------------*/
/* FUNCTIONS ON GAL_T                                    */
/*-------------------------------------------------------*/
gal_t *initGal();
void deallocate_gal(gal_t *thisGal);

/*-------------------------------------------------------*/
/* FUNCTIONS ON GTREE_T                                  */
/*-------------------------------------------------------*/
gtree_t *initGtree(int Ngal);
void reallocGtree(gtree_t **thisGtree, int Ngal);
void deallocate_gtree(gtree_t *thisGtree);
void deallocate_gtreeList(gtree_t **theseGtrees, int numGtrees);

#endif
