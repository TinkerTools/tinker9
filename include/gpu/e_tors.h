#ifndef TINKER_GPU_E_TORS_H_
#define TINKER_GPU_E_TORS_H_


#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
extern int ntors;
extern int (*itors)[4];
extern real (*tors1)[4];
extern real (*tors2)[4];
extern real (*tors3)[4];
extern real (*tors4)[4];
extern real (*tors5)[4];
extern real (*tors6)[4];
extern real torsunit;

extern real* et;
extern real* vir_et;

void etors_data(rc_op op);

void etors(int vers);
TINKER_NAMESPACE_END

#endif
