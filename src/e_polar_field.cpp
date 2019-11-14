#include "e_polar.h"
#include "md.h"


TINKER_NAMESPACE_BEGIN
extern void dfield_ewald_real_acc(real (*)[3], real (*)[3]);
extern void dfield_ewald_real_cu(real (*)[3], real (*)[3]);
void dfield_ewald_real(real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   dfield_ewald_real_cu(field, fieldp);
   return;
#endif
   dfield_ewald_real_acc(field, fieldp);
}


void dfield_ewald(real (*field)[3], real (*fieldp)[3])
{
   device_array::zero(n, field, fieldp);

   dfield_ewald_recip_self(field);
   device_array::copy(n, fieldp, field);

   dfield_ewald_real(field, fieldp);
}


void dfield(real (*field)[3], real (*fieldp)[3])
{
   if (epolar_electyp == elec_t::ewald)
      dfield_ewald(field, fieldp);
   else
      dfield_coulomb(field, fieldp);
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
