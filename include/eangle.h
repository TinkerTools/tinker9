#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"

namespace tinker {
enum class eangle_t : int
{
   in_plane,
   harmonic,
   linear,
   fourier
};

// module angbnd
TINKER_EXTERN int nangle;
TINKER_EXTERN pointer<int, 4> iang;
TINKER_EXTERN pointer<real> ak, anat;

// module angpot
TINKER_EXTERN real angunit;
TINKER_EXTERN real cang, qang, pang, sang;
TINKER_EXTERN pointer<eangle_t> angtyp;

TINKER_EXTERN energy_buffer ea;
TINKER_EXTERN virial_buffer vir_ea;
TINKER_EXTERN grad_prec *deax, *deay, *deaz;
TINKER_EXTERN energy_prec energy_ea;

void eangle_data(rc_op op);

void eangle(int vers);
void eangle_acc(int);
}
