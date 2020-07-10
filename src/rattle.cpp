#include "rattle.h"
#include "mdpq.h"
#include "tool/darray.h"
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
   real dist;
   bool operator<(const HCInfo& h) const
   {
      return this->ir < h.ir;
   }
};


// HC-Molecule is a set of bonds.
using HCMol = std::set<HCInfo>;
bool operator<(const HCMol& h1, const HCMol& h2)
{
   return *h1.begin() < *h2.begin();
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
      nrat = 0;
      darray::deallocate(irat, krat, iratmol);
      darray::deallocate(rattle_xold, rattle_yold, rattle_zold, rattle_moved,
                         rattle_update);
      darray::deallocate(rattle_notdone);
   }


   if (op bitand rc_alloc) {
      rateps = freeze::rateps;
      nrat = freeze::nrat;
      darray::allocate(nrat, &irat, &krat);
      darray::allocate(n, &rattle_xold, &rattle_yold, &rattle_zold,
                       &rattle_moved, &rattle_update);
      darray::allocate(nrat, &rattle_notdone);


      for (int ir = 0; ir < nrat; ++ir) {
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
      nratmol = hc_mols.size();
      darray::allocate(nratmol, &iratmol);


      std::sort(hc_mols.begin(), hc_mols.end());
      std::vector<int> iratm(2 * nratmol);
      std::vector<int> iratn(2 * nrat);
      std::vector<real> kr(nrat);
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


         int mb0 = mbegin;
         for (auto& hc : hcm) {
            iratn[2 * mb0 + 0] = hc.i;
            iratn[2 * mb0 + 1] = hc.k;
            kr[mb0] = hc.dist;
            ++mb0;
         }
      }
      darray::copyin(WAIT_NEW_Q, nrat, irat, iratn.data());
      darray::copyin(WAIT_NEW_Q, nrat, krat, kr.data());
      darray::copyin(WAIT_NEW_Q, nratmol, iratmol, iratm.data());


      hc_mols.clear();
      hc_dict.clear();
   }
}


void rattle(time_prec dt, const pos_prec* xold, const pos_prec* yold,
            const pos_prec* zold)
{
   rattle_acc(dt, xold, yold, zold);
}


void rattle2(time_prec dt, bool do_v)
{
   rattle2_acc(dt, do_v);
}
}
