#include "echglj.h"
#include "box.h"
#include "md.h"
#include "nblist.h"
#include "potent.h"
#include "tool/host_zero.h"
#include "tool/io_fort_str.h"
#include "tool/io_text.h"
#include <array>
#include <map>
#include <set>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/chgpot.hh>
#include <tinker/detail/couple.hh>
#include <tinker/detail/keys.hh>
#include <tinker/detail/kvdws.hh>
#include <tinker/detail/params.hh>
#include <tinker/detail/sizes.hh>
#include <tinker/detail/vdwpot.hh>


namespace tinker {
void echglj_data(rc_op op)
{
   if (!(clist_version() & NBL_SPATIAL))
      return;


   if (!use_potent(charge_term) || !use_potent(vdw_term))
      return;


   if (op & rc_dealloc) {
      ncvexclude = 0;
      vdwpr_in_use = false;
      darray::deallocate(cvexclude, cvexclude_scale);
      darray::deallocate(atom_rad, atom_eps);
   }


   if (op & rc_alloc) {
      // check "VDWPR" keyword
      if (!vdwpr_in_use) {
         auto parse_vdwpr = [](std::string line, int& i, int& k, double& rad,
                               double& eps) {
            try {
               auto vs = Text::split(line);
               std::string ke = vs.at(0);
               Text::upcase(ke);
               if (ke == "VDWPR") {
                  i = std::stoi(vs.at(1));
                  k = std::stoi(vs.at(2));
                  rad = std::stod(vs.at(3));
                  eps = std::stod(vs.at(4));
                  return true;
               }
               return false;
            } catch (...) {
               return false;
            }
         };
         const int* src = atoms::type;
         if (vdwindex == evdw_t::atom_class)
            src = atomid::class_;
         std::set<int> all_tid(src, src + n);
         auto end = all_tid.end();
         int i, k;
         double rad, eps;
         std::string record;
         // first prm
         for (int ii = 0; ii < params::nprm && !vdwpr_in_use; ++ii) {
            fstr_view fsv = params::prmline[ii];
            record = fsv.trim();
            bool okay = parse_vdwpr(record, i, k, rad, eps);
            if (okay) {
               auto iit = all_tid.find(i);
               auto kit = all_tid.find(k);
               if (iit != end && kit != end) {
                  vdwpr_in_use = true;
               }
            }
         }
         // then key
         for (int ii = 0; ii < keys::nkey && !vdwpr_in_use; ++ii) {
            fstr_view fsv = keys::keyline[ii];
            record = fsv.trim();
            bool okay = parse_vdwpr(record, i, k, rad, eps);
            if (okay) {
               auto iit = all_tid.find(i);
               auto kit = all_tid.find(k);
               if (iit != end && kit != end) {
                  vdwpr_in_use = true;
               }
            }
         }
      }


      struct cv
      {
         real c, v;
      };
      auto insert_cv = [](std::map<std::pair<int, int>, cv>& a, int i, int k,
                          real val, char ch) {
         std::pair<int, int> key;
         key.first = i;
         key.second = k;
         auto it = a.find(key);
         if (it == a.end()) {
            cv x;
            x.c = 0;
            x.v = 0;
            if (ch == 'c')
               x.c = val;
            else if (ch == 'v')
               x.v = val;
            a[key] = x;
         } else {
            if (ch == 'c')
               it->second.c = val;
            else if (ch == 'v')
               it->second.v = val;
         }
      };


      // see also attach.h
      const int maxn13 = 3 * sizes::maxval;
      const int maxn14 = 9 * sizes::maxval;
      const int maxn15 = 27 * sizes::maxval;
      const real c2scale = chgpot::c2scale;
      const real c3scale = chgpot::c3scale;
      const real c4scale = chgpot::c4scale;
      const real c5scale = chgpot::c5scale;
      const real v2scale = vdwpot::v2scale;
      const real v3scale = vdwpot::v3scale;
      const real v4scale = vdwpot::v4scale;
      const real v5scale = vdwpot::v5scale;


      std::map<std::pair<int, int>, cv> ik_scale;
      for (int i = 0; i < n; ++i) {
         int nn;
         int bask;


         // c
         if (c2scale != 1) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, c2scale - 1, 'c');
               }
            }
         }
         if (c3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, c3scale - 1, 'c');
               }
            }
         }
         if (c4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, c4scale - 1, 'c');
               }
            }
         }
         if (c5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, c5scale - 1, 'c');
               }
            }
         }


         // v
         if (v2scale != 1 && vdw_exclude_bond == false) {
            nn = couple::n12[i];
            for (int j = 0; j < nn; ++j) {
               int k = couple::i12[i][j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, v2scale - 1, 'v');
               }
            }
         }
         if (v3scale != 1) {
            nn = couple::n13[i];
            bask = i * maxn13;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i13[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, v3scale - 1, 'v');
               }
            }
         }
         if (v4scale != 1) {
            nn = couple::n14[i];
            bask = i * maxn14;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i14[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, v4scale - 1, 'v');
               }
            }
         }
         if (v5scale != 1) {
            nn = couple::n15[i];
            bask = i * maxn15;
            for (int j = 0; j < nn; ++j) {
               int k = couple::i15[bask + j];
               k -= 1;
               if (k > i) {
                  insert_cv(ik_scale, i, k, v5scale - 1, 'v');
               }
            }
         }
      }
      std::vector<int> ik_vec;
      std::vector<real> scal_vec;
      for (auto& it : ik_scale) {
         ik_vec.push_back(it.first.first);
         ik_vec.push_back(it.first.second);
         scal_vec.push_back(it.second.c);
         scal_vec.push_back(it.second.v);
      }
      ncvexclude = ik_scale.size();
      darray::allocate(ncvexclude, &cvexclude, &cvexclude_scale);
      darray::copyin(WAIT_NEW_Q, ncvexclude, cvexclude, ik_vec.data());
      darray::copyin(WAIT_NEW_Q, ncvexclude, cvexclude_scale, scal_vec.data());
      darray::allocate(n, &atom_rad, &atom_eps);
   }


   if (op & rc_init) {
      std::vector<real> vrad(n), veps(n);
      if (vdwindex == evdw_t::atom_type) {
         for (int i = 0; i < n; ++i) {
            int jj = atoms::type[i] - 1;
            vrad[i] = kvdws::rad[jj];
            veps[i] = kvdws::eps[jj];
         }
      } else if (vdwindex == evdw_t::atom_class) {
         for (int i = 0; i < n; ++i) {
            int jj = atomid::class_[i] - 1;
            vrad[i] = kvdws::rad[jj];
            veps[i] = kvdws::eps[jj];
         }
      }
      darray::copyin(WAIT_NEW_Q, n, atom_rad, vrad.data());
      darray::copyin(WAIT_NEW_Q, n, atom_eps, veps.data());
   }
}


void echglj(int vers)
{
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;


   host_zero(energy_ec, virial_ec);
   host_zero(energy_ev, virial_ev);


   assert(vdwtyp == evdw_t::lj);
   assert(radrule == evdw_t::arithmetic);
   assert(epsrule == evdw_t::geometric);
   if (use_ewald()) {
      echarge_ewald_recip_self(vers);
      echglj_rad_arith_eps_geom_ewald_real_cu(vers);
   } else {
      echglj_rad_arith_eps_geom_nonewald_cu(vers);
   }


   if (do_e) {
      if (elrc_vol != 0) {
         energy_prec corr = elrc_vol / volbox();
         energy_ev += corr;
         energy_vdw += corr;
      }
   }
   if (do_v) {
      if (vlrc_vol != 0) {
         virial_prec term = vlrc_vol / volbox();
         virial_ev[0] += term; // xx
         virial_ev[4] += term; // yy
         virial_ev[8] += term; // zz
         virial_vdw[0] += term;
         virial_vdw[4] += term;
         virial_vdw[8] += term;
      }
   }
}
}
