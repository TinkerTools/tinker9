#include "itgiLP22.h"
#include "box.h"
#include "itgEnum.h"
#include "itgiNhc06.h"
#include "lpiston.h"
#include "mathfunc_ou.h"
#include "mdegv.h"
#include "nose.h"
#include "random.h"
#include "tool/io_print.h"
#include "tool/trimatexp.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/mdstuf.hh>
#include <tinker/detail/stodyn.hh>
#include <tinker/detail/units.hh>

namespace tinker {
double* LP22Thermostat::vbarKinetic()
{
   static double ekvbar;
   if (bath::anisotrop) {
      if (box_shape == BoxShape::ORTHO or box_shape == BoxShape::OCT) {
         ekvbar = vbar_matrix[0][0] * vbar_matrix[0][0] +
            vbar_matrix[1][1] * vbar_matrix[1][1] +
            vbar_matrix[2][2] * vbar_matrix[2][2];
      } else if (box_shape == BoxShape::MONO) {
         ekvbar = vbar_matrix[0][0] * vbar_matrix[0][0] +
            vbar_matrix[1][1] * vbar_matrix[1][1] +
            vbar_matrix[2][2] * vbar_matrix[2][2] +
            vbar_matrix[0][2] * vbar_matrix[0][2];
      } else if (box_shape == BoxShape::TRI) {
         ekvbar = vbar_matrix[0][0] * vbar_matrix[0][0] +
            vbar_matrix[1][1] * vbar_matrix[1][1] +
            vbar_matrix[2][2] * vbar_matrix[2][2] +
            vbar_matrix[0][2] * vbar_matrix[0][2] +
            vbar_matrix[0][1] * vbar_matrix[0][1] +
            vbar_matrix[1][2] * vbar_matrix[1][2];
      }
      ekvbar *= 0.5 * qbar;
      return &ekvbar;
   } else {
      return Nhc06Thermostat::vbarKinetic();
   }
}

void LP22Thermostat::scaleVbarVelocity(double scale)
{
   if (bath::anisotrop) {
      if (box_shape == BoxShape::ORTHO or box_shape == BoxShape::OCT) {
         vbar_matrix[0][0] *= scale;
         vbar_matrix[1][1] *= scale;
         vbar_matrix[2][2] *= scale;
      } else if (box_shape == BoxShape::MONO) {
         vbar_matrix[0][0] *= scale;
         vbar_matrix[1][1] *= scale;
         vbar_matrix[2][2] *= scale;
         vbar_matrix[0][2] *= scale;
      } else if (box_shape == BoxShape::TRI) {
         vbar_matrix[0][0] *= scale;
         vbar_matrix[1][1] *= scale;
         vbar_matrix[2][2] *= scale;
         vbar_matrix[0][2] *= scale;
         vbar_matrix[0][1] *= scale;
         vbar_matrix[1][2] *= scale;
      }
   } else {
      Nhc06Thermostat::scaleVbarVelocity(scale);
   }
}

LP22Thermostat::LP22Thermostat()
   : BasicThermostat(ThermostatEnum::m_LP2022)
{
   double dofT = mdstuf::nfree;
   int nhclen = 4;
   int nc = 4;

   // tpart
   m_tpart = new NhcThermostat(nhclen, nc, dofT, NhcThermostat::atomicKinetic,
                               NhcThermostat::scaleAtomicVelocity, "NHC");

   // tbaro
   if (bath::anisotrop) {
      dofT = 1.0;
   } else {
      dofT = 3.0;
      if (box_shape == BoxShape::MONO)
         dofT = 4.0;
      if (box_shape == BoxShape::TRI)
         dofT = 6.0;
   }
   m_tbaro = new NhcThermostat(nhclen, nc, dofT, LP22Thermostat::vbarKinetic,
                               LP22Thermostat::scaleVbarVelocity, "NHCB");
}

LP22Thermostat::~LP22Thermostat()
{
   delete m_tpart;
   delete m_tbaro;
}

void LP22Thermostat::printDetail(FILE* o)
{
   m_tpart->printDetail(o);
   m_tbaro->printDetail(o);
}

void LP22Thermostat::control1b(time_prec dt, bool applyBaro)
{
   if (applyBaro)
      m_tbaro->control1(dt);
   m_tpart->control1(dt);
}

void LP22Thermostat::control2b(time_prec dt, bool applyBaro)
{
   m_tpart->control2(dt, true);
   if (applyBaro)
      m_tbaro->control2(dt, false);
}

//====================================================================//

constexpr int LP22Barostat::m_shapeArray[6][2];

void LP22Barostat::control12Impl(time_prec dt)
{
   double dim = 3.0;
   double kt = units::gasconst * bath::kelvin;
   double b = 2 * m_fric * kt / qbar;
   time_prec dt2 = 0.5 * dt;
   double vol0 = volbox();

   if (bath::anisotrop) {
      double gbar[3][3] = {0};
      for (int i = 0; i < 3; ++i) {
         gbar[i][i] +=
            2 * eksum / m_dofP - vol0 * bath::atmsph / units::prescon;
      }
      for (int k = 0; k < m_arrlen; ++k) {
         int i = m_shapeArray[k][0];
         int j = m_shapeArray[k][1];
         gbar[i][j] += (2 * ekin[i][j] - vir[3 * i + j]);
         gbar[i][j] /= qbar;
         vbar_matrix[i][j] = ornstein_uhlenbeck_process(
            dt2, vbar_matrix[i][j], m_fric, gbar[i][j], b, m_rdn[i][j]);
      }
   } else {
      double al = 1.0 + dim / m_dofP;
      double tr_vir = vir[0] + vir[4] + vir[8];
      double gbar =
         al * 2 * eksum - tr_vir - dim * vol0 * bath::atmsph / units::prescon;
      gbar /= qbar;
      vbar =
         ornstein_uhlenbeck_process(dt2, vbar, m_fric, gbar, b, m_rdn[0][0]);
   }
}

LP22Barostat::LP22Barostat()
   : BasicBarostat(BarostatEnum::LP2022)
{
   m_dofP = mdstuf::nfree;
   m_fric = stodyn::friction;

   m_arrlen = 3;
   if (box_shape == BoxShape::MONO)
      m_arrlen = 4;
   else if (box_shape == BoxShape::TRI)
      m_arrlen = 6;

   double kt = units::gasconst * bath::kelvin;
   qbar = kt * bath::taupres * bath::taupres * m_dofP;
   vbar = 0;
   for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
         vbar_matrix[i][j] = 0;
         m_rdn[i][j] = normal<double>();
      }
   }
}

double LP22Barostat::dof() const
{
   return m_dofP;
}

void LP22Barostat::printDetail(FILE* o)
{
   print(o, " VBar Mass        : %12.2lf\n", qbar);
   printBasic(o);
}

void LP22Barostat::control1(time_prec dt)
{
   control12Impl(dt);
}

void LP22Barostat::control2(time_prec dt)
{
   for (int i = 0; i < 3; ++i) {
      for (int j = i; j < 3; ++j) {
         m_rdn[i][j] = normal<double>();
      }
   }
   control12Impl(dt);
}

void LP22Barostat::control3(time_prec dt)
{
   if (not m_apply)
      return;

   if (bath::anisotrop) {
      double scal[3][3];
      trimat_exp(scal, vbar_matrix, dt);
      double h0[3][3] = {{lvec1.x, lvec1.y, lvec1.z},
                         {lvec2.x, lvec2.y, lvec2.z},
                         {lvec3.x, lvec3.y, lvec3.z}};
      matmul3(h0, scal);
      lvec1.x = h0[0][0], lvec1.y = h0[0][1], lvec1.z = h0[0][2];
      lvec2.x = h0[1][0], lvec2.y = h0[1][1], lvec2.z = h0[1][2];
      lvec3.x = h0[2][0], lvec3.y = h0[2][1], lvec3.z = h0[2][2];
   } else {
      double vt = vbar * dt;
      double eterm2 = std::exp(vt);
      lvec1 *= eterm2;
      lvec2 *= eterm2;
      lvec3 *= eterm2;
   }
   set_default_recip_box();
}

//====================================================================//
}
