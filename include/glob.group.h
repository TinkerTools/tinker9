#pragma once
#include "macro.h"
#include "mdprec.h"

namespace tinker {
/// \brief Partitioning of system into atom groups.
struct Group
{
   /// \brief Total number of atom groups in the system.
   int ngrp;
   /// \brief Contiguous list of the atoms in each group.
   int* kgrp;
   /// \brief Number of the group to which each atom belongs.
   int* grplist;
   /// \brief First and last atom of each group in the list.
   int (*igrp)[2];
   /// \brief Total mass of all the atoms in each group.
   double* grpmass;
   /// \brief Weight for each set of group-group interactions.
   real* wgrp;
};
TINKER_EXTERN Group grp;
}
