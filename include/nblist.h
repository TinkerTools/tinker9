#pragma once
#include "box.h"
#include "dev_array.h"
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
typedef GenericUnit<NBList, GenericUnitVersion::EnableOnDevice> NBListUnit;
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
extern int always_use_nblist;


// size cannot be determined in @c epolar_data(...)
// when neighbor lists are not set up, so these variables
// are declared here
/// index into preconditioner inverse for PCG solver
TINKER_EXTERN device_pointer<int> mindex;
/// preconditioner inverse for induced dipole PCG solver
TINKER_EXTERN device_pointer<real> minv;
TINKER_EXTERN device_pointer<real> minv_exclude_;

void nblist_data(rc_op op);
TINKER_NAMESPACE_END
