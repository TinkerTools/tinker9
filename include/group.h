#pragma once
#include "rc_man.h"


TINKER_NAMESPACE_BEGIN
/// \ingroup md
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
   real* grpmass;
   /// \brief Weight for each set of group-group interactions.
   real* wgrp;


   ~Group();
};
TINKER_EXTERN Group grp;


void group_data(rc_op);
TINKER_NAMESPACE_END
