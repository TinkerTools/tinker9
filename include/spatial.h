#pragma once
#include "gen_unit.h"
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup spatial
 * \brief K-D Tree
 *
 * This is an almost-complete binary tree, where each tree node represents a
 * space partitioning `B`. A node either has two or zero children nodes.
 * A node with no child can be called a leaf. Only one node does not have a
 * parent and that is the root node. For any parent node with its left and
 * right children, regarding their space partitions:
 *    - \f$ B_p = B_l \cup B_r \f$.
 *
 * Every tree is configured with the number of atoms (a variable) and the leaf
 * capacity (a constant). Moreover, the size of an atom block is also set to
 * the leaf capacity, so every leaf node also corresponds to an atom block.
 * Other properties determined by these two parameters are:
 *    - Number of atom blocks;
 *    - Padded number of atoms;
 *    - Number of layers of the tree;
 *    - Number of the tree nodes;
 *    - Etc.
 *
 * The tree nodes can be accessed in two ways:
 *    - By an integer (`node_id`): starting from 0 on the root node, and
 *      iterating the nodes by Breadth First Search (BFS).
 *    - By a pair of integers (`l` and `ml`): level `l` starts from 0 at the
 *      root level; index `ml` restarts from 0 on the first (left most) node of
 *      level `l`; `l` and `ml` can determine a range of atoms in the form of
 *      `[begin, end)`.
 */
struct KDTree
{
   /**
    * \ingroup spatial
    * \brief Space partition.
    */
   struct Bound
   {
      int begin, end;
      real xmin, xmax;
      real ymin, ymax;
      real zmin, zmax;
   };
   // Make sure struct Bound is compact in memory.
   static_assert(sizeof(Bound) == sizeof(int) * 2 + sizeof(real) * 6, "");


   /**
    * \ingroup spatial
    * \brief Atom coordinates for the K-D tree.
    */
   struct Q
   {
      real x, y, z;
   };
   // Make sure struct Q is compact in memory.
   static_assert(sizeof(Q) == sizeof(real) * 3, "");


   /**
    * \brief Leaf capacity and atom block size.
    */
   static constexpr int LEAF = 32;


   /**
    * \brief Number of leaf nodes.
    * \param natoms
    * Number of atoms.
    */
   static int nleaf(int natoms)
   {
      return (natoms + LEAF - 1) / LEAF;
   }


   /**
    * \brief Number of layers.
    * \param natoms
    * Number of atoms.
    * | NLeaves | NLeaves - 1 | NLayers |
    * |:-------:|:-----------:|:-------:|
    * | 1       | 0           | 2       |
    * | 2       | 1           | 2       |
    * | 4       | 3           | 3       |
    * | 5       | 4           | 4       |
    */
   static int nlayer(int natoms)
   {
      return 2 + builtin_floor_log2(nleaf(natoms) - 1);
   }


   /// \brief Number of atoms.
   int n_atom;


   /// \brief Number of atom blocks.
   int n_block;


   /// \brief Padded (by leaf capacity) number of atoms.
   int padded_n_atom;


   /// \brief Number of tree node layers.
   int n_layer;


   Q* xyz;         // of size n
   Bound* bbx;     // of size pow2(n_layer)
   int* reorder;   // of size padded_n
   int* leaf_node; // of size nblock


   ~KDTree();
};
using KDTreeUnit = GenericUnit<KDTree, GenericUnitVersion::EnableOnDevice>;
TINKER_EXTERN KDTreeUnit vtree_unit;
TINKER_EXTERN KDTreeUnit mtree_unit;


void spatial_data(rc_op);
TINKER_NAMESPACE_END
/// \defgroup spatial Spatial Decomposition
