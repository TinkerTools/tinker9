#include "potent.h"
#include <cassert>
#include <tinker/detail/angbnd.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/bndstr.hh>
#include <tinker/detail/improp.hh>
#include <tinker/detail/imptor.hh>
#include <tinker/detail/opbend.hh>
#include <tinker/detail/pitors.hh>
#include <tinker/detail/potent.hh>
#include <tinker/detail/restrn.hh>
#include <tinker/detail/strbnd.hh>
#include <tinker/detail/tors.hh>
#include <tinker/detail/tortor.hh>
#include <tinker/detail/urey.hh>


namespace tinker {
int use_potent(potent_t term)
{
   int val = 0;
   switch (term) {

      // bonded terms

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
   case improp_term:
      val = potent::use_improp;
      break;
   case imptors_term:
      val = potent::use_imptor;
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

      // misc. terms

   case geom_term:
      val = potent::use_geom;
      break;

      // non-bonded terms

   case vdw_term:
      val = potent::use_vdw;
      break;
   case charge_term:
      val = potent::use_charge;
      break;
   case mpole_term:
      val = potent::use_mpole;
      break;
   case polar_term:
      val = potent::use_polar;
      break;

   case chgtrn_term:
      val = potent::use_chgtrn;
      break;
   case disp_term:
      val = potent::use_disp;
      break;
   case repuls_term:
      val = potent::use_repuls;
      break;

   default:
      assert(false);
      break;
   }
   return val;
}


int count_bonded_term(potent_t term)
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
   case improp_term:
      val = improp::niprop;
      break;
   case imptors_term:
      val = imptor::nitors;
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

      // misc. terms

   case geom_term: {
      val = 0;
      val += restrn::npfix;    // position restraints
      val += restrn::ndfix;    // distance restraints
      val += restrn::nafix;    // angle restraints
      val += restrn::ntfix;    // torsional restraints
      val += restrn::ngfix;    // group distance restraints
      val += restrn::nchir;    // chirality restraints
      if (restrn::use_basin) { // Gaussian basin
         val += atoms::n * (atoms::n - 1) / 2;
      }
      if (restrn::use_wall) { // droplet boundary
         val += atoms::n;
      }
   } break;

   default:
      assert(false);
      break;
   }
   return val;
}
}
