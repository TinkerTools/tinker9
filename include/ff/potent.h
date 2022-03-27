#pragma once

namespace tinker {
/// \ingroup ff
enum class Potent
{
   BOND,    // eb
   ANGLE,   // ea
   STRBND,  // eba
   UREY,    // eub
   ANGANG,  // eaa
   OPBEND,  // eopb
   OPDIST,  // eopd
   IMPROP,  // eid
   IMPTORS, // eit
   TORSION, // et
   PITORS,  // ept
   STRTOR,  // ebt
   ANGTOR,  // eat
   TORTOR,  // ett

   VDW,    // ev
   REPULS, // er
   DISP,   // edsp
   CHARGE, // ec
   CHGDPL, // ecd
   DIPOLE, // ed
   MPOLE,  // em
   POLAR,  // ep
   CHGTRN, // ect
   RXNFLD, // erxf
   SOLV,   // es
   METAL,  // elf
   GEOM,   // eg
   EXTRA   // ex
};

/// \ingroup ff
int usePotent(Potent term);
/// \ingroup ff
int countBondedTerm(Potent term);
}
