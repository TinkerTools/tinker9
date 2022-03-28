#pragma once
#include "ff/pchg/echarge.h"
#include "ff/pchg/evdw.h"
#include "tool/rcman.h"

namespace tinker {
void echgljData(RcOp);
void echglj_data_cu(RcOp);
void echglj(int vers);

void echglj_rad_arith_eps_geom_nonewald_cu(int);
void echglj_rad_arith_eps_geom_ewald_real_cu(int);

void pme_stream_start_record_cu(bool use_pmestream);
void pme_stream_start_wait_cu(bool use_pmestream);
void pme_stream_finish_record_cu(bool use_pmestream);
void pme_stream_finish_wait_cu(bool use_pmestream);
#if TINKER_CUDART
inline void pme_stream_start_record(bool use_pmestream)
{
   pme_stream_start_record_cu(use_pmestream);
}
inline void pme_stream_start_wait(bool use_pmestream)
{
   pme_stream_start_wait_cu(use_pmestream);
}
inline void pme_stream_finish_record(bool use_pmestream)
{
   pme_stream_finish_record_cu(use_pmestream);
}
inline void pme_stream_finish_wait(bool use_pmestream)
{
   pme_stream_finish_wait_cu(use_pmestream);
}
#else
inline void pme_stream_start_record(bool use_pmestream) {}
inline void pme_stream_start_wait(bool use_pmestream) {}
inline void pme_stream_finish_record(bool use_pmestream) {}
inline void pme_stream_finish_wait(bool use_pmestream) {}
#endif
}

//====================================================================//
//                                                                    //
//                          Global Variables                          //
//                                                                    //
//====================================================================//

// chglj
namespace tinker {
TINKER_EXTERN int ncvexclude;
TINKER_EXTERN int (*cvexclude)[2];
TINKER_EXTERN real (*cvexclude_scale)[2];

TINKER_EXTERN bool vdwpr_in_use;

TINKER_EXTERN int* mut_coalesced;     // n
TINKER_EXTERN real* chg_coalesced;    // n
TINKER_EXTERN real* radeps_coalesced; // 2*n
}
