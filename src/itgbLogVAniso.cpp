#include "box.h"
#include "integrator.h"
#include "lpiston.h"
#include "mathfunc_ou.h"
#include "nose.h"
#include "random.h"
#include "tool/io.h"
#include "tool/trimatexp.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/units.hh>

namespace tinker {
void AnisoBaroDevice::control_1_2(time_prec dt)
{
   time_prec dt2 = dt * 0.5;
   double vol0 = boxVolume();
   double b = 0;
   if (m_langevin)
      b = std::sqrt(2 * m_fric * units::gasconst * bath::kelvin / qbar);

   double eksu1 = *m_eksum;
   double gbar[3][3] = {0};
   for (int i = 0; i < 3; ++i)
      gbar[i][i] += 2 * eksu1 / dofP - vol0 * bath::atmsph / units::prescon;
   for (int k = 0; k < anisoArrayLength; ++k) {
      int i = anisoArray[k][0];
      int j = anisoArray[k][1];
      gbar[i][j] += (2 * m_ekin[i][j] - m_vir[3 * i + j]);
      gbar[i][j] /= qbar;
      if (m_langevin)
         vbar_matrix[i][j] =
            ornstein_uhlenbeck_process(dt2, vbar_matrix[i][j], m_fric, gbar[i][j], b, m_rdn[i][j]);
      else
         vbar_matrix[i][j] += gbar[i][j] * dt2;
   }
}

AnisoBaroDevice::AnisoBaroDevice(double fric)
   : BasicBarostat()
   , m_vir(nullptr)
   , m_eksum(nullptr)
   , m_ekin(nullptr)
   , m_fric(fric)
   , m_rdn()
   , m_langevin(fric != 0.0)
{
   switch (box_shape) {
   case BoxShape::TRI:
      anisoArrayLength = Tri;
      break;
   case BoxShape::MONO:
      anisoArrayLength = Mono;
      break;
   default:
      anisoArrayLength = OrthoOrOct;
      break;
   }

   if (atomic) {
      m_vir = vir;
      m_eksum = &eksum;
      m_ekin = ekin;
   } else {
      m_vir = lp_vir;
      m_eksum = &lp_eksum;
      m_ekin = lp_ekin;
   }

   if (m_langevin) {
      for (int k = 0; k < anisoArrayLength; ++k) {
         int i = anisoArray[k][0];
         int j = anisoArray[k][1];
         m_rdn[i][j] = normal<double>();
      }
   }

   dofP = mdstuf::nfree;
   double kt = units::gasconst * bath::kelvin;
   qbar = kt * bath::taupres * bath::taupres * dofP;
   vbar = 0;
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         vbar_matrix[i][j] = 0;
      }
   }
}

BarostatEnum AnisoBaroDevice::getBarostatEnum() const
{
   return BarostatEnum::m_LogVAniso;
}

void AnisoBaroDevice::printDetail(FILE* o)
{
   print(o, " VBar Mass          %12.2lf\n", qbar);
   printBasic(o);
}

void AnisoBaroDevice::control1(time_prec dt)
{
   if (not applyBaro)
      return;
   control_1_2(dt);
}

void AnisoBaroDevice::control2(time_prec dt)
{
   if (not applyBaro)
      return;

   for (int k = 0; k < anisoArrayLength; ++k) {
      int i = anisoArray[k][0];
      int j = anisoArray[k][1];
      m_rdn[i][j] = normal<double>();
   }
   control_1_2(dt);
}

void AnisoBaroDevice::control3(time_prec dt)
{
   if (not applyBaro)
      return;

   double scal[3][3];
   trimatExp(scal, vbar_matrix, dt);
   double h0[3][3] = {
      {lvec1.x, lvec1.y, lvec1.z}, {lvec2.x, lvec2.y, lvec2.z}, {lvec3.x, lvec3.y, lvec3.z}};
   matmul3(h0, scal);
   lvec1.x = h0[0][0], lvec1.y = h0[0][1], lvec1.z = h0[0][2];
   lvec2.x = h0[1][0], lvec2.y = h0[1][1], lvec2.z = h0[1][2];
   lvec3.x = h0[2][0], lvec3.y = h0[2][1], lvec3.z = h0[2][2];
   boxSetDefaultRecip();
}
}
