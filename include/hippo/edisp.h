#pragma once
#include "glob/disp.h"
#include "glob/dsppot.h"
#include "tool/rcman.h"

namespace tinker {
bool use_dewald();
void edisp_data(RcOp);
void edisp(int vers);

void edisp_ewald(int vers);
void edisp_ewald_real_cu(int vers);
void edisp_ewald_real_acc(int vers);
void edisp_ewald_recip_self_cu(int vers);
void edisp_ewald_recip_self_acc(int vers);
void disp_pme_conv_acc(int vers);

void edisp_nonewald(int vers);
void edisp_nonewald_cu(int vers);
void edisp_nonewald_acc(int vers);
}
