#ifndef TINKER_GPU_E_VDW_H_
#define TINKER_GPU_E_VDW_H_

#include "util_cxx.h"
#include "util_rc_man.h"

TINKER_NAMESPACE_BEGIN
typedef enum {
  vdw_lj = 0x0001,
  vdw_buck = 0x0002,
  vdw_mm3hb = 0x0004,
  vdw_hal = 0x0008,
  vdw_gauss = 0x0010
} evdw_t;

extern evdw_t vdwtyp;
const char* vdwtyp_str(evdw_t typ);

extern double vdw_switch_cut, vdw_switch_off;

extern real ghal, dhal;
extern real scexp, scalpha;
extern int vcouple;
extern real v2scale, v3scale, v4scale, v5scale;

extern int* ired;
extern real* kred;
extern real *xred, *yred, *zred;

extern int *jvdw, *njvdw;
extern real *radmin, *epsilon;

const int vcouple_decouple = 0;
const int vcouple_annihilate = 1;
extern real* vlam;

extern real* ev;
extern int* nev;
extern real* vir_ev;

void get_evdw_type(evdw_t& typ);
void evdw_data(rc_t rc);

void evdw_reduce_xyz();

void evdw_lj(int vers);
void evdw_buck(int vers);
void evdw_mm3hb(int vers);
void evdw_hal(int vers);
void evdw_gauss(int vers);
void evdw(int vers);
TINKER_NAMESPACE_END

#endif
