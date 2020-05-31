#pragma once
#include "macro.h"
#include "mdprec.h"


namespace tinker {
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
   mass_prec totmass;
   /// \brief Molecular weight for each molecule in the system.
   mass_prec* molmass;
};
TINKER_EXTERN Molecule molecule;
}
