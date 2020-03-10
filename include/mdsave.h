#pragma once
#include "mdprec.h"
#include "rc_man.h"


/**
 * \defgroup md_save  Saving MD Trajectory Snapshots
 * \ingroup md
 */


TINKER_NAMESPACE_BEGIN
void mdsave_async(int istep, time_prec dt);
void mdsave_synchronize();


void mdsave_data(rc_op);
TINKER_NAMESPACE_END
