#ifndef TINKER_E_VDW_H_
#define TINKER_E_VDW_H_

#include "device_vector.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum class evdw_t { lj, buck, mm3hb, hal, gauss };
TINKER_EXTERN evdw_t vdwtyp;

TINKER_EXTERN double vdw_switch_cut, vdw_switch_off;

TINKER_EXTERN real ghal, dhal;
TINKER_EXTERN real scexp, scalpha;
TINKER_EXTERN int vcouple;
TINKER_EXTERN real v2scale, v3scale, v4scale, v5scale;

TINKER_EXTERN DeviceVector<int> ired_vec;
TINKER_EXTERN DeviceVector<real> kred_vec;
TINKER_EXTERN DeviceVector<real> xred_vec, yred_vec, zred_vec;

TINKER_EXTERN DeviceVector<int> jvdw_vec, njvdw_vec;
TINKER_EXTERN DeviceVector<real> radmin_vec, epsilon_vec;

const int vcouple_decouple = 0;
const int vcouple_annihilate = 1;
TINKER_EXTERN DeviceVector<real> vlam_vec;

TINKER_EXTERN int nvdw_excluded_;
TINKER_EXTERN DeviceVector<int, 2> vdw_excluded_vec_;
TINKER_EXTERN DeviceVector<real> vdw_excluded_scale_vec_;

TINKER_EXTERN NonbondedEnergy ev_handle;

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
