#include "ff/amoeba/induce.h"
#include "ff/elec.h"
#include "ff/energy.h"
#include "ff/nblist.h"

namespace tinker {
void dfield_ewald_recip_self(real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_real(real (*field)[3], real (*fieldp)[3]);

void dfield_nonewald_acc(real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_recip_self_acc(real (*field)[3]);
void dfield_ewald_real_acc(real (*field)[3], real (*fieldp)[3]);
void dfield_nonewald_cu(real (*field)[3], real (*fieldp)[3]);
void dfield_ewald_real_cu(real (*field)[3], real (*fieldp)[3]);

void ufield_ewald_recip_self(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);

void ufield_nonewald_acc(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_recip_self_acc(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real_acc(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_nonewald_cu(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
void ufield_ewald_real_cu(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3]);
}

namespace tinker {
void dfield(real (*field)[3], real (*fieldp)[3])
{
   if (useEwald())
      dfieldEwald(field, fieldp);
   else
      dfieldNonEwald(field, fieldp);
}

void dfieldNonEwald(real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      dfield_nonewald_cu(field, fieldp);
   else
#endif
      dfield_nonewald_acc(field, fieldp);
}

void dfieldEwald(real (*field)[3], real (*fieldp)[3])
{
   dfield_ewald_recip_self(field, fieldp);
   dfield_ewald_real(field, fieldp);
}

void dfield_ewald_recip_self(real (*field)[3], real (*fieldp)[3])
{
   dfield_ewald_recip_self_acc(field);
   darray::copy(g::q0, n, fieldp, field);
}

void dfield_ewald_real(real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      dfield_ewald_real_cu(field, fieldp);
   else
#endif
      dfield_ewald_real_acc(field, fieldp);
}

void ufield(const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
   if (useEwald())
      ufieldEwald(uind, uinp, field, fieldp);
   else
      ufieldNonEwald(uind, uinp, field, fieldp);
}

void ufieldNonEwald(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      ufield_nonewald_cu(uind, uinp, field, fieldp);
   else
#endif
      ufield_nonewald_acc(uind, uinp, field, fieldp);
}

void ufieldEwald(const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
   ufield_ewald_recip_self(uind, uinp, field, fieldp);
   ufield_ewald_real(uind, uinp, field, fieldp);
}

void ufield_ewald_recip_self(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
   ufield_ewald_recip_self_acc(uind, uinp, field, fieldp);
}

void ufield_ewald_real(
   const real (*uind)[3], const real (*uinp)[3], real (*field)[3], real (*fieldp)[3])
{
#if TINKER_CUDART
   if (mlist_version() & NBL_SPATIAL)
      ufield_ewald_real_cu(uind, uinp, field, fieldp);
   else
#endif
      ufield_ewald_real_acc(uind, uinp, field, fieldp);
}
}
