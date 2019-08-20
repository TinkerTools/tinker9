#ifndef TINKER_E_VDW_H_
#define TINKER_E_VDW_H_

#include "list.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum class evdw_t { lj, buck, mm3hb, hal, gauss };
TINKER_EXTERN evdw_t vdwtyp;

TINKER_EXTERN double vdw_switch_cut, vdw_switch_off;

TINKER_EXTERN real ghal, dhal;
TINKER_EXTERN real scexp, scalpha;
TINKER_EXTERN int vcouple;
TINKER_EXTERN real v2scale, v3scale, v4scale, v5scale;

TINKER_EXTERN int* ired;
TINKER_EXTERN real* kred;
TINKER_EXTERN real *xred, *yred, *zred;

TINKER_EXTERN int *jvdw, *njvdw;
TINKER_EXTERN real *radmin, *epsilon;

const int vcouple_decouple = 0;
const int vcouple_annihilate = 1;
TINKER_EXTERN real* vlam;

TINKER_EXTERN int* nev;
TINKER_EXTERN real_buffer_t* ev;
TINKER_EXTERN real_buffer_t* vir_ev;
TINKER_EXTERN size_t bufsize_ev;

void evdw_data(rc_op op);

void evdw_reduce_xyz();

void evdw_lj(int vers);
void evdw_buck(int vers);
void evdw_mm3hb(int vers);
void evdw_hal(int vers);
void evdw_gauss(int vers);
void evdw(int vers);
TINKER_NAMESPACE_END

#endif
