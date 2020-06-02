#pragma once
#include "mdprec.h"
#include "rc_man.h"


namespace tinker {
void mdsave_async(int istep, time_prec dt);
void mdsave_synchronize();


void mdsave_data(rc_op);
}
