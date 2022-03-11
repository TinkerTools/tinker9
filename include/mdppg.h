#pragma once
#include "mdprec.h"

namespace tinker {
void velAvbfIso(int nrespa, vel_prec a, vel_prec b, const grad_prec* gx1, const grad_prec* gy1,
                const grad_prec* gz1, const grad_prec* gx2, const grad_prec* gy2,
                const grad_prec* gz2);
void velAvbfAni(int nrespa, vel_prec a[3][3], vel_prec b[3][3], const grad_prec* gx1,
                const grad_prec* gy1, const grad_prec* gz1, const grad_prec* gx2,
                const grad_prec* gy2, const grad_prec* gz2);
}
