#include "integrator.h"
#include "lpiston.h"
#include "md.h"
#include "nose.h"
#include <tinker/detail/mdstuf.hh>

namespace tinker {
Nhc06Thermostat::Nhc06Thermostat()
   : BasicThermostat()
{
   double dofT;
   int nhclen = 4;
   int nc = 4;

   // tpart
   if (atomic) {
      dofT = mdstuf::nfree;
      m_tpart = new NhcDevice(
         nhclen, nc, dofT, NhcDevice::kineticAtomic, NhcDevice::scaleVelocityAtomic, "NHC");
   } else {
      dofT = 3.0 * (rattle_dmol.nmol - 1);
      m_tpart = new NhcDevice(nhclen, nc, dofT, Nhc06Thermostat::kineticRattleGroup,
         Nhc06Thermostat::scaleVelocityRattleGroup, "NHC");
   }

   // tbaro
   m_tbaro = new NhcDevice(nhclen, nc, dofVbar(), Nhc06Thermostat::kineticVbar,
      Nhc06Thermostat::scaleVelocityVbar, "NHCBaro");
}

Nhc06Thermostat::~Nhc06Thermostat()
{
   delete m_tbaro;
   delete m_tpart;
}

void Nhc06Thermostat::printDetail(FILE* o)
{
   m_tpart->printDetail(o);
   m_tbaro->printDetail(o);
}

void Nhc06Thermostat::control1(time_prec dt)
{
   if (applyBaro)
      m_tbaro->control1(dt);
   m_tpart->control1(dt);
}

void Nhc06Thermostat::control2(time_prec dt, bool save)
{
   m_tpart->control2(dt, save);
   if (applyBaro)
      m_tbaro->control2(dt, false);
}

double Nhc06Thermostat::kineticRattleGroup()
{
   lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
   lp_mol_kinetic();
   return lp_eksum;
}

void Nhc06Thermostat::scaleVelocityRattleGroup(double scale)
{
   double s2 = scale * scale;
   lp_eksum *= s2;
   for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
         lp_ekin[i][j] *= s2;
   lp_propagate_mol_vel(scale - 1);
}

double Nhc06Thermostat::kineticVbar()
{
   double ekvbar;
   int i, j, k;
   ekvbar = 0;
   switch (anisoArrayLength) {
   case Tri:
   case Mono:
   case OrthoOrOct:
   case SemiIso:
      for (k = 0; k < anisoArrayLength; ++k) {
         i = anisoArray[k][0];
         j = anisoArray[k][1];
         ekvbar += 0.5 * qbar * vbar_matrix[i][j];
      }
      break;
   default:
      ekvbar += 0.5 * qbar * vbar * vbar;
      break;
   }
   return ekvbar;
}

void Nhc06Thermostat::scaleVelocityVbar(double scale)
{
   int i, j, k;
   switch (anisoArrayLength) {
   case Tri:
   case Mono:
   case OrthoOrOct:
      for (k = 0; k < anisoArrayLength; ++k) {
         i = anisoArray[k][0];
         j = anisoArray[k][1];
         vbar_matrix[i][j] *= scale;
      }
      break;
   case SemiIso:
      for (k = 0; k < OrthoOrOct; ++k) {
         i = anisoArray[k][0];
         j = anisoArray[k][1];
         vbar_matrix[i][j] *= scale;
      }
      break;
   default:
      vbar *= scale;
      break;
   }
}

double Nhc06Thermostat::dofVbar()
{
   double dof;
   switch (anisoArrayLength) {
   case Tri:
      dof = static_cast<double>(Tri);
      break;
   case Mono:
      dof = static_cast<double>(Mono);
      break;
   case OrthoOrOct:
      dof = static_cast<double>(OrthoOrOct);
      break;
   case SemiIso:
      dof = static_cast<double>(SemiIso);
      break;
   default:
      dof = 1.0;
      break;
   }
   return dof;
}
}
