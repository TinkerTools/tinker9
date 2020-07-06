#pragma once
#include "mod.mutant.h"
#include "mod.vdw.h"
#include "mod.vdwpot.h"
#include "tool/rc_man.h"


/**
 * \todo Test lj, buck, mm3hb, gauss.
 */


namespace tinker {
void evdw_data(rc_op op);
/**
 * \ingroup vdw
 * \brief Lennard-Jones 12-6 potential.
 *
 * \f[ U(r|r_m,\epsilon) = \epsilon [(r_m/r)^{12} - 2(r_m/r)^6] \f]
 */
void elj(int vers);
void elj_acc(int);
void elj_cu(int);
void elj14_cu(int);
void ebuck(int vers);
void ebuck_acc(int);
void emm3hb(int vers);
void emm3hb_acc(int);
/**
 * \ingroup vdw
 * \brief
 * [Halgren buffered 14-7 potential.](https://doi.org/10.1021/ja00046a032)
 *
 * \f[
 * U(r|r_m,\epsilon) = \epsilon
 * \left(\frac{1+\gamma}{\rho^m+\gamma}-2\right)
 * \left(\frac{1+\delta}{\rho+\delta}\right)^{n-m}
 * \f]
 * \f[ \rho = r/r_m,\ n=14,\ m=7,\ \gamma=0.12,\ \delta=0.07 \f]
 *
 * [Soft core buffered 14-7 potential.](https://doi.org/10.1002/jcc.21681)
 *
 * \f[
 * U(r|\lambda) = \lambda^t\epsilon
 * \left(\frac{1.12}{\alpha(1-\lambda)^2+\rho^7+0.12}-2\right)
 * \frac{1.07^5}{\alpha(1-\lambda)^2+(\rho+0.07)^7}
 * \f]
 * \f[ t=5,\ \alpha=0.7 \f]
 */
void ehal(int vers);
void ehal_acc(int);
void ehal_cu(int);
void ehal_reduce_xyz();
void ehal_resolve_gradient();
void ehal_reduce_xyz_acc();
void ehal_resolve_gradient_acc();
void egauss(int vers);
void egauss_acc(int);
void evdw(int vers);
}
