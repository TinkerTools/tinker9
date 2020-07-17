#pragma once
#include "glob.cudalib.h"
#include "tool/rc_man.h"


namespace tinker {
void cudalib_data(rc_op);
void stream2_begin();
// synchronize events in other streams with the `nonblk` stream
void stream2_sync();
}
