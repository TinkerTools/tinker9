#pragma once
#include "mod.disp.h"
#include "mod.dsppot.h"
#include "tool/rc_man.h"


namespace tinker {
bool use_dewald();
void edisp_data(rc_op);
void edisp(int vers);
void edisp_cu(int vers);
}
