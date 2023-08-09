#ifndef COPY_HALO_GAL_H
#define COPY_HALO_GAL_H

void copy_halo_to_gal(halo_t *thisHalos, gal_t *thisGal);
gal_t *gal_from_halo(halo_t *thisHalo);
gtree_t *gtree_from_tree(tree_t *thisTree);
int32_t gtrees_from_trees(tree_t ***thisTreeList, int numTree, gtree_t ***thisGtreeList);

#endif
