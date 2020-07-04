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
}
