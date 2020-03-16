#ifndef HALO_TREE_H
#define HALO_TREE_H

typedef struct{
  int64_t ID;
  int64_t descID;
  int32_t localID;
  int32_t localDescID;
  int32_t numProg;

  int32_t snapnumber;
  float scalefactor;
  float descScalefactor;
  
  float pos[3];
  float vel[3];
  float Mvir;
  float Rvir;
  float velDisp;
  float velMax;
  float halfmassRadius;
  float spin;
  float spinBullock;
  float scalefactorLastMajorMerger;
}halo_t;

typedef struct{
  int numHalos;
  halo_t *halos;
}tree_t;

/*-------------------------------------------------------*/
/* FUNCTIONS ON HALO_T                                   */
/*-------------------------------------------------------*/
halo_t *initHalo();
void deallocate_halo(halo_t *thisOutHalo);

/*-------------------------------------------------------*/
/* FUNCTIONS ON TREE_T                                   */
/*-------------------------------------------------------*/
tree_t *initTree(int Nhalos);
void reallocTree(tree_t **thisTree, int Nhalos);
void deallocate_tree(tree_t *thisTree);
void deallocate_treeList(tree_t **theseTrees, int numTrees);

#endif