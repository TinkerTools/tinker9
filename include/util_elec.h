#ifndef TINKER_UTIL_ELEC_H_
#define TINKER_UTIL_ELEC_H_

#include "mod_elec.h"
#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
int use_elec();
void elec_data(rc_t rc);
void elec_init(int vers);
void torque(int vers);
TINKER_NAMESPACE_END

#endif
