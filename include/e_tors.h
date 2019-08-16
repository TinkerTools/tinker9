#ifndef TINKER_E_TORS_H_
#define TINKER_E_TORS_H_

#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
TINKER_EXTERN int ntors;
TINKER_EXTERN int (*itors)[4];
TINKER_EXTERN real (*tors1)[4];
TINKER_EXTERN real (*tors2)[4];
TINKER_EXTERN real (*tors3)[4];
TINKER_EXTERN real (*tors4)[4];
TINKER_EXTERN real (*tors5)[4];
TINKER_EXTERN real (*tors6)[4];
TINKER_EXTERN real torsunit;

TINKER_EXTERN real* et;
TINKER_EXTERN real* vir_et;

void etors_data(rc_op op);

void etors(int vers);
TINKER_NAMESPACE_END

#endif
