#include "e_polar.h"
#include "md.h"
#include <tinker/detail/bound.hh>


TINKER_NAMESPACE_BEGIN
void dfield_ewald_real(real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   extern void dfield_ewald_real_cu(real(*)[3], real(*)[3]);
   dfield_ewald_real_cu(field, fieldp);
#else
   extern void dfield_ewald_real_acc(real(*)[3], real(*)[3]);
   dfield_ewald_real_acc(field, fieldp);
#endif
}


void dfield_ewald(real (*field)[3], real (*fieldp)[3])
{
   device_array::zero(DM_NEW_Q_PROCEED, n, field, fieldp);
   dfield_ewald_recip_self(field);
   device_array::copy(false, n, fieldp, field);
   dfield_ewald_real(field, fieldp);
}


void dfield_coulomb(real (*field)[3], real (*fieldp)[3])
{
   extern void dfield_coulomb_acc(real(*)[3], real(*)[3]);
#if TINKER_CUDART
   extern void dfield_coulomb_cu(real(*)[3], real(*)[3]);
   if (bound::use_bounds)
      dfield_coulomb_cu(field, fieldp);
   else
#endif
      dfield_coulomb_acc(field, fieldp);
}


void dfield(real (*field)[3], real (*fieldp)[3])
{
   if (epolar_electyp == elec_t::ewald)
      dfield_ewald(field, fieldp);
   else
      dfield_coulomb(field, fieldp);
}


void ufield_ewald_real(const real (*uind)[3], const real (*uinp)[3],
                       real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   extern void ufield_ewald_real_cu(const real(*)[3], const real(*)[3],
                                    real(*)[3], real(*)[3]);
   ufield_ewald_real_cu(uind, uinp, field, fieldp);
#else
   extern void ufield_ewald_real_acc(const real(*)[3], const real(*)[3],
                                     real(*)[3], real(*)[3]);
   ufield_ewald_real_acc(uind, uinp, field, fieldp);
#endif
}


void ufield_ewald(const real (*uind)[3], const real (*uinp)[3],
                  real (*field)[3], real (*fieldp)[3])
{
   device_array::zero(DM_NEW_Q_PROCEED, n, field, fieldp);
   ufield_ewald_recip_self(uind, uinp, field, fieldp);
   ufield_ewald_real(uind, uinp, field, fieldp);
}


void ufield_coulomb(const real (*uind)[3], const real (*uinp)[3],
                    real (*field)[3], real (*fieldp)[3])
{
   extern void ufield_coulomb_acc(const real(*)[3], const real(*)[3],
                                  real(*)[3], real(*)[3]);
#if TINKER_CUDART
   extern void ufield_coulomb_cu(const real(*)[3], const real(*)[3], real(*)[3],
                                 real(*)[3]);
   if (bound::use_bounds)
      ufield_coulomb_cu(uind, uinp, field, fieldp);
   else
#endif
      ufield_coulomb_acc(uind, uinp, field, fieldp);
}


void ufield(const real (*uind)[3], const real (*uinp)[3], real (*field)[3],
            real (*fieldp)[3])
{
   if (epolar_electyp == elec_t::ewald)
      ufield_ewald(uind, uinp, field, fieldp);
   else
      ufield_coulomb(uind, uinp, field, fieldp);
}
TINKER_NAMESPACE_END
