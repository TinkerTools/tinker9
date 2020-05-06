#pragma once
#include "darray.h"
#include "energy_buffer.h"
#include "rc_man.h"


/**
 * \todo Test lj, buck, mm3hb, gauss.
 */


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup vdw
 * \brief Constant flags for the VDW energy functions.
 */
enum class evdw_t
{
   lj,    ///< Lennard-Jones 12-6 potential.
   buck,  ///< Buckingham potential.
   mm3hb, ///< MM3 exp-6 potential.
   hal,   ///< Halgren buffered 14-7 potential.
   gauss, ///< Gaussian expansion VDW potential.


   decouple = 0,   ///< VDW lambda type: decouple.
   annihilate = 1, ///< VDW lambda type: annihilate.
};
TINKER_EXTERN evdw_t vdwtyp;
TINKER_EXTERN evdw_t vcouple;


/**
 * \ingroup vdw
 * \brief Value of \f$ \gamma \f$ in buffered 14-7 vdw potential.
 */
TINKER_EXTERN real ghal;
/**
 * \ingroup vdw
 * \brief Value of \f$ \delta \f$ in buffered 14-7 vdw potential.
 */
TINKER_EXTERN real dhal;
/**
 * \ingroup vdw
 * \brief Exponential factor for soft core buffered 14-7 potential.
 */
TINKER_EXTERN real scexp;
/**
 * \ingroup vdw
 * \brief Scale factor \f$ \alpha \f$ for soft core buffered 14-7 potential.
 */
TINKER_EXTERN real scalpha;


TINKER_EXTERN real v2scale, v3scale, v4scale, v5scale;
extern bool vdw_exclude_bond;


/**
 * \ingroup vdw
 * \brief Halgren buffered 14-7 reduced x, y, z coordinates for each atom.
 */
TINKER_EXTERN pointer<real> xred, yred, zred;
/**
 * \ingroup vdw
 * \brief Halgren buffered 14-7 reduced vdw gradients for each atom.
 */
TINKER_EXTERN pointer<grad_prec> gxred, gyred, gzred;
TINKER_EXTERN pointer<int> ired;
TINKER_EXTERN pointer<real> kred;


/**
 * \ingroup vdw
 * \brief Number of unique values in the #jvdw array.
 */
TINKER_EXTERN int njvdw;
/**
 * \ingroup vdw
 * \brief Type or class index into vdw parameters for each atom.
 * The indices have been sorted and start from 0.
 */
TINKER_EXTERN pointer<int> jvdw;


/**
 * \ingroup vdw
 * \brief Minimum energy distance (#radmin) or well depth parameter (#epsilon)
 * for each #jvdw pair. Element `[j1][j2]` is accessed by `[njvdw*j1 + j2]`.
 * \see njvdw
 */
TINKER_EXTERN pointer<real> radmin, epsilon;


/**
 * \ingroup vdw
 * \brief VDW 1-4 parameters: minimum energy distance and well depth.
 * \see radmin epsilon
 */
TINKER_EXTERN pointer<real> radmin4, epsilon4;
TINKER_EXTERN int nvdw14;
TINKER_EXTERN pointer<int, 2> vdw14ik;


/**
 * \ingroup vdw
 * \brief
 * State weighting values (lambda) of all atoms for van der Waals potentials.
 */
TINKER_EXTERN pointer<int> mut;
TINKER_EXTERN real vlam;


TINKER_EXTERN int nvexclude;
TINKER_EXTERN pointer<int, 2> vexclude;
TINKER_EXTERN pointer<real> vexclude_scale;


TINKER_EXTERN count_buffer nev;
TINKER_EXTERN energy_buffer ev;
TINKER_EXTERN virial_buffer vir_ev;
TINKER_EXTERN grad_prec *devx, *devy, *devz;
TINKER_EXTERN energy_prec energy_ev;


/**
 * \ingroup vdw
 * \brief Long-range energy correction (lrc), used as `e += lrc/volume`.
 * \note Must be 0 if system is unbound.
 */
TINKER_EXTERN energy_prec elrc_vol;
/**
 * \ingroup vdw
 * \brief Long-range virial correction (lrc), used as `v(i,i) += lrc/volume`.
 * \note Must be 0 if system is unbound.
 */
TINKER_EXTERN virial_prec vlrc_vol;


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
TINKER_NAMESPACE_END
