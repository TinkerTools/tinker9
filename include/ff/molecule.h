#pragma once
#include "ff/precision.h"
#include "tool/rcman.h"

namespace tinker {
/// \addtogroup ff
/// \{

/// Individual molecules in current system.
struct Molecule
{
   int nmol;        ///< Total number of separate molecules in the system.
   int (*imol)[2];  ///< First and last atom of each molecule in the list.
   int* kmol;       ///< Contiguous list of the atoms in each molecule.
   int* molecule;   ///< Number of the molecule to which each atom belongs.
   double totmass;  ///< Total weight of all the molecules in the system.
   double* molmass; ///< Molecular weight for each molecule in the system.
};

/// Partitioning of system into atom groups.
struct Group
{
   int ngrp;        ///< Total number of atom groups in the system.
   int* kgrp;       ///< Contiguous list of the atoms in each group.
   int* grplist;    ///< Number of the group to which each atom belongs.
   int (*igrp)[2];  ///< First and last atom of each group in the list.
   double* grpmass; ///< Total mass of all the atoms in each group.
   real* wgrp;      ///< Weight for each set of group-group interactions.
};

void coupleData(RcOp);
void moleculeData(RcOp);
void groupData(RcOp);

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

const int couple_maxn12 = 8;
TINKER_EXTERN int (*couple_i12)[couple_maxn12];
TINKER_EXTERN int* couple_n12;
TINKER_EXTERN Molecule molecule;
TINKER_EXTERN Group grp;

/// \}
}
