#include "field_aplus.h"
#include "elec.h"
#include "field.h"
#include "md.h"
#include "nblist.h"


namespace tinker {
void dfield_aplus(real (*field)[3])
{
   if (use_ewald())
      dfield_aplus_ewald(field);
   else
      dfield_aplus_nonewald(field);
}


void dfield_aplus_nonewald(real (*field)[3])
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      dfield_aplus_nonewald_cu(field);
   else
#endif
      dfield_aplus_nonewald_acc(field);
}


void dfield_aplus_ewald(real (*field)[3])
{
   dfield_aplus_ewald_recip_self(field);
   dfield_aplus_ewald_real(field);
}


void dfield_aplus_ewald_recip_self(real (*field)[3])
{
   dfield_ewald_recip_self_acc(field);
}


void dfield_aplus_ewald_real(real (*field)[3])
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      dfield_aplus_ewald_real_cu(field);
   else
#endif
      dfield_aplus_ewald_real_acc(field);
}


void ufield_aplus(const real (*uind)[3], real (*field)[3])
{
   if (use_ewald())
      ufield_aplus_ewald(uind, field);
   else
      ufield_aplus_nonewald(uind, field);
}


void ufield_aplus_nonewald(const real (*uind)[3], real (*field)[3])
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      ufield_aplus_nonewald_cu(uind, field);
   else
#endif
      ufield_aplus_nonewald_acc(uind, field);
}


void ufield_aplus_ewald(const real (*uind)[3], real (*field)[3])
{
   ufield_aplus_ewald_recip_self(uind, field);
   ufield_aplus_ewald_real(uind, field);
}


void ufield_aplus_ewald_recip_self(const real (*uind)[3], real (*field)[3])
{
   ufield_aplus_ewald_recip_self_acc(uind, field);
}


void ufield_aplus_ewald_real(const real (*uind)[3], real (*field)[3])
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      ufield_aplus_ewald_real_cu(uind, field);
   else
#endif
      ufield_aplus_ewald_real_acc(uind, field);
}
}
