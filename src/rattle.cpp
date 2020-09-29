#include "rattle.h"
#include "mdegv.h"
#include "mdpq.h"
#include "tool/darray.h"
#include "tool/energy_buffer.h"
#include "tool/io_print.h"
#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <tinker/detail/freeze.hh>


namespace tinker {
bool use_rattle()
{
   return freeze::use_rattle;
}


namespace {
// holonomic constraint information
struct HCInfo
{
   int i, k; // atoms i and k
   int ir;   // value from 0 to nrat-1
   pos_prec dist;
   bool operator<(const HCInfo& h) const
   {
      return this->ir < h.ir;
   }


   // sorted by (i,k) pair
   struct less_ver2
   {
      bool operator()(const HCInfo& h1, const HCInfo& h2) const
      {
         if (h1.i < h2.i) {
            return true;
         } else if (h1.i > h2.i) {
            return false;
         } else {
            return h1.k < h2.k;
         }
      }
   };
};


// HC-Molecule is a set of bonds.
using HCMol = std::set<HCInfo>;
bool operator<(const HCMol& h1, const HCMol& h2)
{
   return *h1.begin() < *h2.begin();
}


bool HCMol_is_water(const HCMol& h, int& a, int& b, int& c, pos_prec& ab,
                    pos_prec& ac, pos_prec& bc)
{
   if (h.size() != 3)
      return false;
   // h has 3 constraints
   std::vector<HCInfo> m(h.begin(), h.end());
   std::sort(m.begin(), m.end(), HCInfo::less_ver2());
   if (m[0].i == m[1].i and m[0].k == m[2].i and m[1].k == m[2].k) {
      a = m[0].i;
      b = m[0].k;
      c = m[1].k;
      ab = m[0].dist;
      ac = m[1].dist;
      bc = m[2].dist;
      return true;
   } else {
      return false;
   }
}


// all of the HC-Molecules
std::vector<HCMol> hc_mols;
// hc_dict[i] gives where to find atom i
std::map<int, size_t> hc_dict;
}


void rattle_data(rc_op op)
{
   if (not use_rattle())
      return;


   if (op bitand rc_dealloc) {
      rateps = 0;
      nratwt = 0;
      darray::deallocate(iratwt, kratwt);
      nrat = 0;
      nratmol = 0;
      darray::deallocate(irat, krat, iratmol);
      darray::deallocate(rattle_xold, rattle_yold, rattle_zold);
      // darray::deallocate(rattle_moved, rattle_update);
      // darray::deallocate(rattle_bigdelta);
   }


   if (op bitand rc_alloc) {
      // save data from Fortran library
      for (int ir = 0; ir < freeze::nrat; ++ir) {
         int i0 = freeze::irat[ir * 2 + 0] - 1;
         int k0 = freeze::irat[ir * 2 + 1] - 1;
         HCInfo hc;
         hc.i = std::min(i0, k0);
         hc.k = std::max(i0, k0);
         hc.ir = ir;
         hc.dist = freeze::krat[ir];


         // check if atoms i or k are already recorded
         auto iloc = hc_dict.find(hc.i);
         auto kloc = hc_dict.find(hc.k);
         if (iloc == hc_dict.end() && kloc == hc_dict.end()) {
            // both i and k are new
            HCMol hcm;
            hcm.insert(hc);         // add this HC to a new HC-Molecule
            hc_mols.push_back(hcm); // record this new HC-Molecule
            size_t sz = hc_mols.size() - 1;
            hc_dict[hc.i] = sz; // record atom i
            hc_dict[hc.k] = sz; // record atom k
         } else if (iloc == hc_dict.end()) {
            // k is recorded, i is new
            size_t sz = kloc->second;
            hc_dict[hc.i] = sz; // record atom i
            auto& hcm = hc_mols[sz];
            hcm.insert(hc);
         } else if (kloc == hc_dict.end()) {
            // i is recorded, k is new
            size_t sz = iloc->second;
            hc_dict[hc.k] = sz; // record atom k
            auto& hcm = hc_mols[sz];
            hcm.insert(hc);
         } else {
            assert(iloc->second == kloc->second);
            size_t sz = iloc->second;
            auto& hcm = hc_mols[sz];
            hcm.insert(hc);
         }
      }


      rateps = freeze::rateps;


      // find water-like constraints in hc_mols
      std::vector<int> veciwater;
      std::vector<pos_prec> veckwater;
      for (auto& it : hc_mols) {
         int a, b, c;
         pos_prec ab, ac, bc;
         if (HCMol_is_water(it, a, b, c, ab, ac, bc)) {
            veciwater.push_back(a);
            veciwater.push_back(b);
            veciwater.push_back(c);
            veckwater.push_back(ab);
            veckwater.push_back(ac);
            veckwater.push_back(bc);
            it.clear();
         }
      }
      assert(veciwater.size() % 3 == 0);
      nratwt = veciwater.size() / 3;
      darray::allocate(nratwt, &iratwt, &kratwt);
      darray::copyin(PROCEED_NEW_Q, nratwt, iratwt, veciwater.data());
      darray::copyin(WAIT_NEW_Q, nratwt, kratwt, veckwater.data());


      // erase water-like constraints in hc_mols
      hc_mols.erase(
         std::remove_if(hc_mols.begin(), hc_mols.end(),
                        [](const HCMol& h) { return h.size() == 0; }),
         hc_mols.end());
      std::sort(hc_mols.begin(), hc_mols.end());


      nratmol = hc_mols.size();
      std::vector<int> iratm(2 * nratmol);
      std::vector<int> iratn;
      std::vector<pos_prec> kr;
      for (int i = 0; i < nratmol; ++i) {
         const auto& hcm = hc_mols[i];
         int msize = hcm.size();
         int mbegin;
         if (i == 0)
            mbegin = 0;
         else
            mbegin = iratm[2 * i - 1];
         iratm[2 * i + 0] = mbegin;
         iratm[2 * i + 1] = mbegin + msize;


         for (auto& hc : hcm) {
            iratn.push_back(hc.i);
            iratn.push_back(hc.k);
            kr.push_back(hc.dist);
         }
      }
      nrat = kr.size();


      darray::allocate(nrat, &irat, &krat);
      darray::allocate(nratmol, &iratmol);
      darray::copyin(WAIT_NEW_Q, nrat, irat, iratn.data());
      darray::copyin(WAIT_NEW_Q, nrat, krat, kr.data());
      darray::copyin(WAIT_NEW_Q, nratmol, iratmol, iratm.data());


      darray::allocate(n, &rattle_xold, &rattle_yold, &rattle_zold);
      // darray::allocate(n, &rattle_moved, &rattle_update);
      // darray::allocate(nrat, &rattle_bigdelta);


      hc_mols.clear();
      hc_dict.clear();
   }
}


void rattle(time_prec dt, const pos_prec* xold, const pos_prec* yold,
            const pos_prec* zold)
{
   rattle_settle_acc(dt, xold, yold, zold);
   rattle_acc(dt, xold, yold, zold);
}


void rattle2(time_prec dt, bool do_v)
{
   if (do_v) {
      darray::zero(PROCEED_NEW_Q, buffer_size(), vir_buf);
   }


   rattle2_settle_acc(dt, do_v);
   rattle2_acc(dt, do_v);


   if (do_v) {
      virial_prec v[9];
      virial_reduce(v, vir_buf);
      for (int iv = 0; iv < 9; ++iv) {
         vir[iv] += v[iv];
      }
   }
}

// Langevin Piston
void rattle_lp(time_prec dt,
              pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
              const pos_prec* xold, const pos_prec* yold,const pos_prec* zold)
{
   rattle_settle_acc_lp(dt, xnew, ynew, znew, xold, yold, zold);
   rattle_acc_lp(dt, xnew, ynew, znew, xold, yold, zold);
}

void rattle2_lf(time_prec dt,
                pos_prec* vx_lp, pos_prec* vy_lp, pos_prec* vz_lp,
                const vel_prec* vx_new, const vel_prec* vy_new, const vel_prec* zx_new,
                const pos_prec* xold, const pos_prec* yold, const pos_prec* zold)
{
   darray::zero(PROCEED_NEW_Q, buffer_size(), vir_buf);

   rattle2_lf_acc(dt, vx_lp, vy_lp, vz_lp, vx_new, vy_new, zx_new, xold, yold, zold);

   virial_prec v[9];
   virial_reduce(v, vir_buf);
   for (int iv = 0; iv < 9; ++iv) {
      vir[iv] += v[iv];
   }
}

}
