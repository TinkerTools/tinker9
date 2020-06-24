#pragma once
#include "mod.disp.h"
#include "mod.dsppot.h"
#include "tool/rc_man.h"


namespace tinker {
bool use_dewald();
void edisp_data(rc_op);
void edisp(int vers);


void edisp_ewald(int vers);
void edisp_ewald_real_cu(int vers);
void edisp_ewald_recip_self_cu(int vers);
void disp_pme_conv_acc(int vers);


void edisp_nonewald(int vers);
void edisp_nonewald_cu(int vers);
}
