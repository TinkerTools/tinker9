#pragma once
#include "ff/precision.h"
#include "tool/rcman.h"

namespace tinker {
/// \ingroup ff
/// \brief Individual molecules in current system.
struct Molecule
{
   /// \brief Total number of separate molecules in the system.
   int nmol;
   /// \brief First and last atom of each molecule in the list.
   int (*imol)[2];
   /// \brief Contiguous list of the atoms in each molecule.
   int* kmol;
   /// \brief Number of the molecule to which each atom belongs.
   int* molecule;
   /// \brief Total weight of all the molecules in the system.
   double totmass;
   /// \brief Molecular weight for each molecule in the system.
   double* molmass;
};

/// \ingroup ff
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
}

namespace tinker {
/// \ingroup ff
void coupleData(RcOp);
/// \ingroup ff
void moleculeData(RcOp);
/// \ingroup ff
void groupData(RcOp);
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

namespace tinker {
/// \ingroup ff
const int couple_maxn12 = 8;
/// \ingroup ff
TINKER_EXTERN int (*couple_i12)[couple_maxn12];
/// \ingroup ff
TINKER_EXTERN int* couple_n12;

/// \ingroup ff
TINKER_EXTERN Molecule molecule;

/// \ingroup ff
TINKER_EXTERN Group grp;
}
