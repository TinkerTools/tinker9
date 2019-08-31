#ifndef TINKER_E_ANGLE_H_
#define TINKER_E_ANGLE_H_

#include "device_vector.h"
#include "energy_buffer.h"
#include "rc_man.h"

TINKER_NAMESPACE_BEGIN
enum class eangle_t : int { in_plane, harmonic, linear, fourier };

// module angbnd
TINKER_EXTERN int nangle;
TINKER_EXTERN DeviceVector<int, 4> iang_vec;
TINKER_EXTERN DeviceVector<real> ak_vec, anat_vec;

// module angpot
TINKER_EXTERN real angunit;
TINKER_EXTERN real cang, qang, pang, sang;
TINKER_EXTERN DeviceVector<eangle_t> angtyp_vec;

TINKER_EXTERN BondedEnergy ea_handle;

void eangle_data(rc_op op);

void eangle(int vers);
TINKER_NAMESPACE_END

#endif
