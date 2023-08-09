#ifndef OUTGAL_H
#define OUTGAL_H

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
#ifndef FIRST
  float fracMgasMer;
#endif
#if defined WITHMETALS
  float MgasNew;
  float MgasEj;
#endif
  float Mgas;
  float Mstar;
  
#ifndef FIRST
  float fesc;
  float Nion;
  float fej;
#endif
  float feff;
  float fg;
  float zreion;
  float photHI_bg;
  
#if defined WITHMETALS
  float Mmetal[3];
  float MmetalIni[3];
  float MmetalNew[3];
  float fracMmetalMer[3];
  float MmetalEj[3];
  
  float Mdust;
#endif
}outgal_t; 

typedef struct{
  int32_t numGal;
  outgal_t *galaxies;
}outgtree_t;

/*-------------------------------------------------------*/
/* FUNCTIONS ON OUTGAL_T                                 */
/*-------------------------------------------------------*/

outgal_t *initOutGalTree();
void deallocate_outgaltree(outgal_t *thisOutGal);

/*-------------------------------------------------------*/
/* FUNCTIONS ON GTREE_T                                  */
/*-------------------------------------------------------*/

outgtree_t *initOutGtree(int numGal);
void reallocOutGtree(outgtree_t **thisOutGtree, int numGal);
void deallocate_outgtree(outgtree_t *thisOutGtree);
void deallocate_outgtreeList(outgtree_t **theseOutGtrees, int numOutGtrees);

#endif
