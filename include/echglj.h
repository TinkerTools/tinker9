#pragma once
#include "echarge.h"
#include "evdw.h"
#include "glob.accasync.h"
#include "glob.chglj.h"
#include "tool/rc_man.h"


namespace tinker {
void echglj_data(rc_op op);
void echglj_data_cu(rc_op);
void echglj_cu_sync_pme_stream(bool use_pmestream);
void echglj(int vers);


void echglj_rad_arith_eps_geom_nonewald_cu(int);
void echglj_rad_arith_eps_geom_ewald_real_cu(int);


extern real* chg_coalesced;    // n
extern real* radeps_coalesced; // 2*n
}
