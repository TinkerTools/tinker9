#include "potent.h"
#include <cassert>
#include <ext/tinker/detail/angbnd.hh>
#include <ext/tinker/detail/bndstr.hh>
#include <ext/tinker/detail/opbend.hh>
#include <ext/tinker/detail/pitors.hh>
#include <ext/tinker/detail/potent.hh>
#include <ext/tinker/detail/strbnd.hh>
#include <ext/tinker/detail/tors.hh>
#include <ext/tinker/detail/tortor.hh>
#include <ext/tinker/detail/urey.hh>

TINKER_NAMESPACE_BEGIN
int use_potent (potent_t term)
{
   int val = 0;
   switch (term) {

      // bonded term

   case bond_term:
      val = potent::use_bond;
      break;
   case angle_term:
      val = potent::use_angle;
      break;
   case strbnd_term:
      val = potent::use_strbnd;
      break;
   case urey_term:
      val = potent::use_urey;
      break;
   case opbend_term:
      val = potent::use_opbend;
      break;
   case torsion_term:
      val = potent::use_tors;
      break;
   case pitors_term:
      val = potent::use_pitors;
      break;
   case tortor_term:
      val = potent::use_tortor;
      break;

      // non-bonded term

   case vdw_term:
      val = potent::use_vdw;
      break;
   case mpole_term:
      val = potent::use_mpole;
      break;
   case polar_term:
      val = potent::use_polar;
      break;
   default:
      assert (false);
      break;
   }
   return val;
}

int count_bonded_term (potent_t term)
{
   int val = -1;
   switch (term) {
   case bond_term:
      val = bndstr::nbond;
      break;
   case angle_term:
      val = angbnd::nangle;
      break;
   case strbnd_term:
      val = strbnd::nstrbnd;
      break;
   case urey_term:
      val = urey::nurey;
      break;
   case opbend_term:
      val = opbend::nopbend;
      break;
   case torsion_term:
      val = tors::ntors;
      break;
   case pitors_term:
      val = pitors::npitors;
      break;
   case tortor_term:
      val = tortor::ntortor;
      break;
   default:
      assert (false);
      break;
   }
   return val;
}
TINKER_NAMESPACE_END
