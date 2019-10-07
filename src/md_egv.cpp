
#include "md.h"
#include <map>

TINKER_NAMESPACE_BEGIN
static bool use_ev_ ()
{
   return rc_flag & (calc::analyz | calc::energy | calc::virial);
}

static void grad_data_ (rc_op op)
{
   if (!(rc_flag & calc::grad))
      return;

   if (op & rc_dealloc) {
      device_array::deallocate (gx, gy, gz);
   }

   if (op & rc_alloc) {
      device_array::allocate (n, &gx, &gy, &gz);
   }

   // we can never assume whether or not deriv::desum was allocated, because it
   // was allocated inside subroutine gradient(...), which would be skipped in
   // subroutine mdinit() if a dyn file existed to restart a simulation.

   // if (op & rc_init) {
   // copy in deriv::sum to gx, gy, and gz
   // }
}

void zero_egv ()
{
   zero_egv (rc_flag & calc::vmask);
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
static void ev_data_ (rc_op op)
{
   if (!use_ev_ ())
      return;

   if (op & rc_dealloc) {
      Count::clear ();
      Energy::clear ();
      Virial::clear ();

      esum_handle.close ();
      vir_handle.close ();
   }

   if (op & rc_alloc) {
      if (!(rc_flag & calc::analyz)) {
         esum_handle = Energy::open ();
         vir_handle = Virial::open ();
      }
   }
}

void egv_data (rc_op op)
{
   rc_man ev42_{ev_data_, op};
   rc_man grad42_{grad_data_, op};
}
TINKER_NAMESPACE_END

TINKER_NAMESPACE_BEGIN
void zero_egv (int vers)
{
   if (vers & calc::analyz) {
      for (int i = 0; i < Count::size (); ++i) {
         Count u = i;
         u->zero ();
      }
   }

   if (vers & calc::energy) {
      for (int i = 0; i < Energy::size (); ++i) {
         Energy u = i;
         u->zero ();
      }
   }

   if (vers & calc::virial) {
      for (int i = 0; i < Virial::size (); ++i) {
         Virial u = i;
         u->zero ();
      }
   }

   if (vers & calc::grad) {
      device_array::zero (n, gx, gy, gz);
   }
}

void sum_energies (int vers)
{
   if (vers & calc::energy) {
      if (rc_flag & calc::analyz) {
         esum = 0;
         for (int i = 0; i < Energy::size (); ++i) {
            Energy u = i;
            real e;
            u->sum (&e);
            esum += e;
         }
      } else {
         esum_handle->sum (&esum);
      }
   }

   if (vers & calc::virial) {
      if (rc_flag & calc::analyz) {
         for (int iv = 0; iv < 9; ++iv)
            vir[iv] = 0;
         for (int i = 0; i < Virial::size (); ++i) {
            Virial u = i;
            real v[9];
            u->sum (v);
            for (int iv = 0; iv < 9; ++iv)
               vir[iv] += v[iv];
         }
      } else {
         vir_handle->sum (&vir[0]);
      }
   }
}
TINKER_NAMESPACE_END
