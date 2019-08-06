#ifndef TINKER_NBLIST_H_
#define TINKER_NBLIST_H_

#include "gen_unit.h"
#include "rc_man.h"
#include "rt.h"

TINKER_NAMESPACE_BEGIN
/// @brief
/// pairwise neighbor list indices and storage
struct NBList {
  enum {
    double_loop = 1, ///< double loop version
    nblist = 2       ///< neighbor list version
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

typedef GenericUnit<NBList, 1> NBListUnit;
TINKER_EXTERN NBListUnit vlist_unit;
TINKER_EXTERN NBListUnit dlist_unit;
TINKER_EXTERN NBListUnit clist_unit;
TINKER_EXTERN NBListUnit mlist_unit;
TINKER_EXTERN NBListUnit ulist_unit;

void nblist_data(rc_op op);
TINKER_NAMESPACE_END

#endif
