#include "energy.h"
#include "io_fort_str.h"
#include "md.h"
#include "nblist.h"
#include "platform.h"
#include "spatial.h"
#include <cassert>
#include <tinker/detail/bath.hh>
#include <tinker/detail/inform.hh>
#include <tinker/detail/mdstuf.hh>

TINKER_NAMESPACE_BEGIN
void md_data(rc_op op)
{
   if ((calc::md & rc_flag) == 0)
      return;

   rc_man intg42_{integrate_data, op};
   rc_man save42_{mdsave_data, op};
}

static void (*intg)(int, real);

void integrate_data(rc_op op)
{
   if (op & rc_dealloc) {
      if (intg == respa_fast_slow)
         device_array::deallocate(gx1, gy1, gz1, gx2, gy2, gz2);

      intg = nullptr;
   }

   if (op & rc_init) {
      if (bath::isothermal) {
         fstr_view th = bath::thermostat;
         if (th == "BERENDSEN")
            thermostat = thermo_berendsen;
         else if (th == "BUSSI")
            thermostat = thermo_bussi;
         else if (th == "ANDERSEN")
            thermostat = thermo_andersen;
         else if (th == "NOSE-HOOVER")
            thermostat = thermo_nose_hoover_chain;
         else
            assert(false);
      } else {
         thermostat = thermo_null;
      }

      if (bath::isobaric) {
         fstr_view br = bath::barostat;
         if (br == "BERENDSEN")
            barostat = baro_berendsen;
         else if (br == "BUSSI")
            barostat = baro_bussi;
         else if (br == "NOSE-HOOVER")
            barostat = baro_nose_hoover_chain;
         else if (br == "MONTECARLO")
            barostat = baro_montecarlo;
         else
            assert(false);
      } else {
         barostat = baro_null;
      }

      fstr_view itg = mdstuf::integrate;
      intg = nullptr;
      if (itg == "VERLET") {
         intg = velocity_verlet;
         // need full gradient to start/restart the simulation
         energy(rc_flag);
      } else if (itg == "STOCHASTIC") {
      } else if (itg == "BAOAB") {
      } else if (itg == "BUSSI") {
      } else if (itg == "NOSE-HOOVER") {
      } else if (itg == "GHMC") {
      } else if (itg == "RIGIDBODY") {
      } else if (itg == "RESPA") {
         intg = respa_fast_slow;
         // need fast and slow gradients to start/restart the simulation
         device_array::allocate(n, &gx1, &gy1, &gz1, &gx2, &gy2, &gz2);
         // save fast gradients to gx1 etc.
         energy(rc_flag, RESPA_FAST, respa_tsconfig());
         copy_energy(rc_flag, nullptr, gx1, gy1, gz1, nullptr);
         // save slow gradients to gx2 etc.
         energy(rc_flag, RESPA_SLOW, respa_tsconfig());
         copy_energy(rc_flag, nullptr, gx2, gy2, gz2, nullptr);
      } else {
         // beeman
         assert(false);
      }
   }
}

void kinetic(real& temp)
{
   extern void kinetic_acc(real&);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM) {
      extern void kinetic_cu(real&);
      kinetic_cu(temp);
   } else
#endif
      kinetic_acc(temp);
}

extern void thermo_bussi_acc(real dt, real temp);
void temper(real dt, real& temp)
{
   kinetic(temp);
   if (thermostat == thermo_null)
      return;

   if (thermostat == thermo_bussi)
      thermo_bussi_acc(dt, temp);
   else
      assert(false);
}

extern void mdrest_acc(int istep);
void mdrest(int istep)
{
   mdrest_acc(istep);
}

void propagate_xyz(real dt, int check_nblist)
{
   extern void propagate_xyz_acc(real);
   propagate_xyz_acc(dt);
   if (check_nblist)
      nblist_data(rc_evolve);
}

void propagate_velocity(real dt, const real* grx, const real* gry,
                        const real* grz)
{
   extern void propagate_velocity_acc(real, const real*, const real*,
                                      const real*);
   propagate_velocity_acc(dt, grx, gry, grz);
}

void propagate_velocity2(real dt, const real* grx, const real* gry,
                         const real* grz, real dt2, const real* grx2,
                         const real* gry2, const real* grz2)
{
   extern void propagate_velocity2_acc(real, const real*, const real*,
                                       const real*, real, const real*,
                                       const real*, const real*);
   propagate_velocity2_acc(dt, grx, gry, grz, dt2, grx2, gry2, grz2);
}

void propagate(int nsteps, real dt_ps)
{
   for (int istep = 1; istep <= nsteps; ++istep) {
      intg(istep, dt_ps);

      // mdstat
      if (istep % inform::iwrite == 0) {
         real temp;
         kinetic(temp);
         mdsave_async(istep, dt_ps);
      }
      mdrest(istep);
   }
   mdsave_synchronize();
}
TINKER_NAMESPACE_END
