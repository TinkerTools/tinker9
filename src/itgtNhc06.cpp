#include "itgtNhc06.h"
#include "lpiston.h"
#include "mdpq.h"
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
      m_tpart =
         new NhcThermostat(nhclen, nc, dofT, NhcThermostat::kineticAtomic,
                           NhcThermostat::scaleVelocityAtomic, "NHC");
   } else {
      dofT = 3.0 * (rattle_dmol.nmol - 1);
      m_tpart = new NhcThermostat(
         nhclen, nc, dofT, Nhc06Thermostat::kineticRattleGroup,
         Nhc06Thermostat::scaleVelocityRattleGroup, "NHC");
   }

   // tbaro
   m_tbaro =
      new NhcThermostat(nhclen, nc, dofVbar(), Nhc06Thermostat::kineticVbar,
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

double* Nhc06Thermostat::kineticRattleGroup()
{
   lp_center_of_mass(vx, vy, vz, ratcom_vx, ratcom_vy, ratcom_vz);
   lp_mol_kinetic();
   return &lp_eksum;
}

void Nhc06Thermostat::scaleVelocityRattleGroup(double scale)
{
   __PlaceHolderMessage("Impl pending...");
}

double* Nhc06Thermostat::kineticVbar()
{
   static double ekvbar;
   int i, j, k;
   ekvbar = 0;
   switch (anisoArrayLength) {
   case Tri:
   case Mono:
   case OrthoOrOct:
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
   return &ekvbar;
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
      dof = 6.0;
      break;
   case Mono:
      dof = 4.0;
      break;
   case OrthoOrOct:
      dof = 3.0;
      break;
   default:
      dof = 1.0;
      break;
   }
   return dof;
}
}
