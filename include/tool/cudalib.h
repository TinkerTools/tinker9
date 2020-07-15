#pragma once
#include "glob.cudalib.h"
#include "tool/rc_man.h"


namespace tinker {
void cudalib_data(rc_op);
// synchronize events in other streams with the `nonblk` stream
void sync_events_with_nonblk();
}
