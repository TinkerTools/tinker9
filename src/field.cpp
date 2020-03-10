#include "field.h"
#include "elec.h"
#include "md.h"
#include "nblist.h"


TINKER_NAMESPACE_BEGIN
void dfield(real (*field)[3], real (*fieldp)[3])
{
   if (use_ewald())
      dfield_ewald(field, fieldp);
   else
      dfield_nonewald(field, fieldp);
}


void dfield_nonewald(real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlist_version() == NBL_SPATIAL)
      dfield_nonewald_cu(field, fieldp);
   else
#endif
      dfield_nonewald_acc(field, fieldp);
}


void dfield_ewald(real (*field)[3], real (*fieldp)[3])
{
   dfield_ewald_recip_self(field, fieldp);
   dfield_ewald_real(field, fieldp);
}


void dfield_ewald_recip_self(real (*field)[3], real (*fieldp)[3])
{
   dfield_ewald_recip_self_acc(field);
   darray::copy(PROCEED_NEW_Q, n, fieldp, field);
}


void dfield_ewald_real(real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlist_version() == NBL_SPATIAL)
      dfield_ewald_real_cu(field, fieldp);
   else
#endif
      dfield_ewald_real_acc(field, fieldp);
}


void ufield(const real (*uind)[3], const real (*uinp)[3], real (*field)[3],
            real (*fieldp)[3])
{
   if (use_ewald())
      ufield_ewald(uind, uinp, field, fieldp);
   else
      ufield_nonewald(uind, uinp, field, fieldp);
}


void ufield_nonewald(const real (*uind)[3], const real (*uinp)[3],
                     real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlist_version() == NBL_SPATIAL)
      ufield_nonewald_cu(uind, uinp, field, fieldp);
   else
#endif
      ufield_nonewald_acc(uind, uinp, field, fieldp);
}


void ufield_ewald(const real (*uind)[3], const real (*uinp)[3],
                  real (*field)[3], real (*fieldp)[3])
{
   ufield_ewald_recip_self(uind, uinp, field, fieldp);
   ufield_ewald_real(uind, uinp, field, fieldp);
}


void ufield_ewald_recip_self(const real (*uind)[3], const real (*uinp)[3],
                             real (*field)[3], real (*fieldp)[3])
{
   ufield_ewald_recip_self_acc(uind, uinp, field, fieldp);
}


void ufield_ewald_real(const real (*uind)[3], const real (*uinp)[3],
                       real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlist_version() == NBL_SPATIAL)
      ufield_ewald_real_cu(uind, uinp, field, fieldp);
   else
#endif
      ufield_ewald_real_acc(uind, uinp, field, fieldp);
}
TINKER_NAMESPACE_END
