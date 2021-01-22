#pragma once
#include "macro.h"
#include "tool/energy_buffer.h"


namespace tinker {
TINKER_EXTERN int npfix;
TINKER_EXTERN int* ipfix;
TINKER_EXTERN int (*kpfix)[3];
TINKER_EXTERN real* xpfix;
TINKER_EXTERN real* ypfix;
TINKER_EXTERN real* zpfix;
TINKER_EXTERN real (*pfix)[2];


/// \ingroup geom
/// \brief Number of group distance restraints to be applied.
TINKER_EXTERN int ngfix;
/// \ingroup geom
/// \brief Group numbers defining each group distance restraint.
TINKER_EXTERN int (*igfix)[2];
/// \ingroup geom
/// \brief Force constant and target range for each group distance.
TINKER_EXTERN real (*gfix)[3];


TINKER_EXTERN int ndfix;
TINKER_EXTERN int (*idfix)[2];
TINKER_EXTERN real (*dfix)[3];


TINKER_EXTERN int nafix;
TINKER_EXTERN int (*iafix)[3];
TINKER_EXTERN real (*afix)[3];


TINKER_EXTERN int ntfix;
TINKER_EXTERN int (*itfix)[4];
TINKER_EXTERN real (*tfix)[3];


TINKER_EXTERN energy_buffer eg;
TINKER_EXTERN virial_buffer vir_eg;
TINKER_EXTERN grad_prec* degx;
TINKER_EXTERN grad_prec* degy;
TINKER_EXTERN grad_prec* degz;
TINKER_EXTERN energy_prec energy_eg;
TINKER_EXTERN virial_prec virial_eg[9];
}
