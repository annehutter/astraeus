#ifndef GAL_GTREE_H
#define GAL_GTREE_H

typedef struct{
#if defined WITHROCKSTARID
  int64_t ID;
  int64_t descID;
#endif
 
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

  float MgasIni;
  float fracMgasMer;
#if defined WITHMETALS
  float MgasNew;
  float MgasEj;
#endif
  float Mgas;
  float Mstar;

  float fesc;
  float Nion;
  float fej;
  float feff;
  float fg;
  float photHI_bg;
  float zreion;
  
  float *stellarmasshistory;

#if defined WITHMETALS
  float Mmetal[3];
  float MmetalIni[3];
  float MmetalNew[3];
  float fracMmetalMer[3];
  float MmetalEj[3];
  float igmMetallicity[3];

  float *metalmasshistory;
  
  float Mdust;
  float MdustEj;
  float igmDustFraction;
#endif
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
