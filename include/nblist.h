#pragma once
#include "box.h"
#include "darray.h"
#include "gen_unit.h"
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
/// \brief
/// Verlet list: pairwise neighbor list indices and storage.
struct NBList
{
   enum
   {
      double_loop = 0x1, ///< Double loop version.
      nblist = 0x2,      ///< Neighbor list version.
      spatial = 0x4      ///< Spatial decomposition version.
   };

   int* nlst;     ///< number of sites in list for each atom
   int* lst;      ///< all of the sites in list
   int* update;   ///< update flag for each atom
   real* xold;    ///< old x coordinates
   real* yold;    ///< old y coordinates
   real* zold;    ///< old z coordinates
   const real* x; ///< current x coordinates
   const real* y; ///< current y coordinates
   const real* z; ///< current z coordinates
   int maxnlst;   ///< max number of neighbors for each atom
   real cutoff;   ///< list cutoff distance
   real buffer;   ///< width of the neighbor list buffer region

   ~NBList();
};
using NBListUnit = GenericUnit<NBList, GenericUnitVersion::EnableOnDevice>;
TINKER_EXTERN NBListUnit vlist_unit;
TINKER_EXTERN NBListUnit dlist_unit;
TINKER_EXTERN NBListUnit clist_unit;
TINKER_EXTERN NBListUnit mlist_unit;
TINKER_EXTERN NBListUnit ulist_unit;


int vlist_version();
int dlist_version();
int clist_version();
int mlist_version();
int ulist_version();


void nblist_data(rc_op op);


void nblist_build_acc(NBListUnit); // rc_init
void nblist_update_acc(NBListUnit);


void refresh_neighbors();
TINKER_NAMESPACE_END
