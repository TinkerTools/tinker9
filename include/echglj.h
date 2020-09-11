#pragma once
#include "echarge.h"
#include "evdw.h"
#include "glob.chglj.h"
#include "tool/rc_man.h"


namespace tinker {
void echglj_data(rc_op op);
void echglj(int vers);


void echglj_rad_arith_eps_geom_nonewald_cu(int);
void echglj_rad_arith_eps_geom_ewald_real_cu(int);


extern real* chg_coalesced;    // n
extern real* radeps_coalesced; // 2*n
// #define TINKER_ECHGLJ_USE_COALESCED_GRAD 0
#define TINKER_ECHGLJ_USE_COALESCED_GRAD 1
extern grad_prec* gx_coalesced;
extern grad_prec* gy_coalesced;
extern grad_prec* gz_coalesced;
}
