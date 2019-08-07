#ifndef TINKER_E_BOND_H_
#define TINKER_E_BOND_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum { bond_harmonic = 0x001, bond_morse = 0x002 };
TINKER_EXTERN int bndtyp;

TINKER_EXTERN real cbnd, qbnd, bndunit;
TINKER_EXTERN int nbond;
TINKER_EXTERN int (*ibnd)[2];
TINKER_EXTERN real *bl, *bk;

TINKER_EXTERN real* eb;
TINKER_EXTERN real* vir_eb;

void ebond_data(rc_op op);

void ebond(int vers);
TINKER_NAMESPACE_END

#endif
