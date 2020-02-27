#pragma once
#include "dev_array.h"
#include "energy_buffer.h"
#include "rc_man.h"

/**
 * \todo Test lj, buck, mm3hb, gauss, and mutant.
 * \todo Add vdw correction.
 */

TINKER_NAMESPACE_BEGIN
/**
 * \ingroup vdw
 * \brief Constant flags for the VDW energy functions.
 */
enum class evdw_t
{
   lj,
   buck,
   mm3hb,
   /// Halgren buffered 14-7 potential.
   hal,
   gauss,

   decouple = 0,
   annihilate = 1,
};
TINKER_EXTERN evdw_t vdwtyp;

/// \ingroup vdw
/// Value of \f$ \gamma \f$ in buffered 14-7 vdw potential.
TINKER_EXTERN real ghal;
/// \ingroup vdw
/// Value of \f$ \delta \f$ in buffered 14-7 vdw potential.
TINKER_EXTERN real dhal;
/// \ingroup vdw
/// Exponential factor for soft core buffered 14-7 potential.
TINKER_EXTERN real scexp;
/// \ingroup vdw
/// Scale factor \f$ \alpha \f$ for soft core buffered 14-7 potential.
TINKER_EXTERN real scalpha;
/// \ingroup vdw
/// Van der Waals lambda type.
/// \see evdw_t::decouple
/// \see evdw_t::annihilate
TINKER_EXTERN evdw_t vcouple;
TINKER_EXTERN real v2scale, v3scale, v4scale, v5scale;

TINKER_EXTERN device_pointer<int> ired;
TINKER_EXTERN device_pointer<real> kred;

/// \ingroup vdw
/// Reduced x, y, z coordinates for each atom and the vdw gradients on each
/// reduced site.
/// \{
TINKER_EXTERN device_pointer<real> xred, yred, zred;
TINKER_EXTERN device_pointer<real> gxred, gyred, gzred;
/// \}

/**
 * \ingroup vdw
 * Number of unique values in the \c jvdw array.
 * \see jvdw
 */
TINKER_EXTERN int njvdw;
/** \ingroup vdw
 * Type or class index into vdw parameters for each atom.
 * The indices have been sorted and start from 0.
 */
TINKER_EXTERN device_pointer<int> jvdw;
/// \ingroup vdw
/// \{
/// Minimum energy distance or well depth parameter for each \c jvdw pair.
/// Element `[j1][j2]` is accessed by `[njvdw*j1 + j2]`.
/// \see njvdw
TINKER_EXTERN device_pointer<real> radmin, epsilon;
/// \}

/// \ingroup vdw
/// State weighting values \f$ \lambda \f$ of all atoms for van der Waals
/// potentials.
TINKER_EXTERN device_pointer<real> vlam;

TINKER_EXTERN int nvexclude_;
TINKER_EXTERN device_pointer<int, 2> vexclude_;
TINKER_EXTERN device_pointer<real> vexclude_scale_;

TINKER_EXTERN count_buffer nev;
TINKER_EXTERN energy_buffer ev;
TINKER_EXTERN virial_buffer vir_ev;

/// \ingroup vdw
/// \brief Long-range energy correction (lrc), used as `e += lrc/pbc_volume`.
TINKER_EXTERN real elrc_vol;
/// \ingroup vdw
/// \brief
/// Long-range virial correction (lrc), used as `v(i,i) += lrc/pbc_volume`.
TINKER_EXTERN real vlrc_vol;

void evdw_data(rc_op op);

void evdw_reduce_xyz();
void evdw_resolve_gradient();

void evdw_lj(int vers);
void evdw_buck(int vers);
void evdw_mm3hb(int vers);
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
void evdw_hal(int vers);
void evdw_gauss(int vers);
void evdw(int vers);
TINKER_NAMESPACE_END
