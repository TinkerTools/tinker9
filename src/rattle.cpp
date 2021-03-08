#include "rattle.h"
#include "mdcalc.h"
#include "mdegv.h"
#include "mdpq.h"
#include "tool/darray.h"
#include "tool/energy_buffer.h"
#include "tool/io_print.h"
#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <tinker/detail/atomid.hh>
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


namespace {
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

// methine group, -S-H, etc.
bool HCMol_is_ch(const HCMol& h, int& a, int& b, pos_prec& ab)
{
   if (h.size() != 1)
      return false;
   std::vector<HCInfo> m(h.begin(), h.end());
   a = m[0].i;
   b = m[0].k;
   ab = m[0].dist;
   return true;
}


// methylene group, -NH2, etc.
bool HCMol_is_ch2(const HCMol& h, int& a, int& b, int& c, pos_prec& ab,
                  pos_prec& ac)
{
   if (h.size() != 2)
      return false;
   std::vector<HCInfo> m(h.begin(), h.end());
   std::sort(m.begin(), m.end(), HCInfo::less_ver2());
   auto& m0 = m[0];
   auto& m1 = m[1];
   if (m0.i == m1.i) {
      a = m0.i;
      b = m0.k;
      c = m1.k;
      ab = m0.dist;
      ac = m1.dist;
      return true;
   } else if (m0.i == m1.k) {
      a = m0.i;
      b = m0.k;
      c = m1.i;
      ab = m0.dist;
      ac = m1.dist;
      return true;
   } else if (m0.k == m1.i) {
      a = m0.k;
      b = m0.i;
      c = m1.k;
      ab = m0.dist;
      ac = m1.dist;
      return true;
   } else if (m0.k == m1.k) {
      a = m0.k;
      b = m0.i;
      c = m1.i;
      ab = m0.dist;
      ac = m1.dist;
      return true;
   }
   return false;
}


// methyl group
bool HCMol_is_ch3(const HCMol& h, int& a, int& b, int& c, int& d, pos_prec& ab,
                  pos_prec& ac, pos_prec& ad)
{
   if (h.size() != 3)
      return false;
   std::vector<HCInfo> m(h.begin(), h.end());
   std::sort(m.begin(), m.end(), HCInfo::less_ver2());
   auto& m0 = m[0];
   auto& m1 = m[1];
   auto& m2 = m[2];
   if (((m0.i == m1.i) or (m0.i == m1.k)) and
       ((m0.i == m2.i) or (m0.i == m2.k))) {
      a = m0.i;
      b = m0.k;
      if (m0.i == m1.i)
         c = m1.k;
      else
         c = m1.i;
      if (m0.i == m2.i)
         d = m2.k;
      else
         d = m2.i;
      ab = m0.dist;
      ac = m1.dist;
      ad = m2.dist;
      return true;
   } else if (((m0.k == m1.i) or (m0.k == m1.k)) and
              ((m0.k == m2.i) or (m0.k == m2.k))) {
      a = m0.k;
      b = m0.i;
      if (m0.k == m1.i)
         c = m1.k;
      else
         c = m1.i;
      if (m0.k == m2.i)
         d = m2.k;
      else
         d = m2.i;
      ab = m0.dist;
      ac = m1.dist;
      ad = m2.dist;
      return true;
   }
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


   if (op & rc_dealloc) {
      rattle_dmol.nmol = 0;
      rattle_dmol.totmass = 0;
      darray::deallocate(rattle_dmol.imol, rattle_dmol.kmol,
                         rattle_dmol.molecule, rattle_dmol.molmass);
      if (rc_flag & calc::md) {
         darray::deallocate(ratcom_x, ratcom_y, ratcom_z, ratcom_px, ratcom_py,
                            ratcom_pz, ratcom_gx, ratcom_gy, ratcom_gz,
                            ratcom_massfrac);
         ratcom_x = nullptr;
         ratcom_y = nullptr;
         ratcom_z = nullptr;
         ratcom_px = nullptr;
         ratcom_py = nullptr;
         ratcom_pz = nullptr;
         ratcom_gx = nullptr;
         ratcom_gy = nullptr;
         ratcom_gz = nullptr;
         ratcom_massfrac = nullptr;
      }

      rateps = 0;
      nratwt = 0;
      darray::deallocate(iratwt, kratwt);
      nratch = 0;
      darray::deallocate(iratch, kratch);
      nratch2 = 0;
      darray::deallocate(iratch2, kratch2);
      nratch3 = 0;
      darray::deallocate(iratch3, kratch3);
      nrat = 0;
      nratmol = 0;
      darray::deallocate(irat, krat, iratmol);
      darray::deallocate(rattle_xold, rattle_yold, rattle_zold);
   }


   if (op & rc_alloc) {
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


      // find the "rattle-molecules"
      {
         // host vectors
         std::vector<int> hvec_imol, hvec_kmol, hvec_molec;
         std::vector<mass_prec> hvec_molmass, hvec_massfrac;

         // first fill hvec_kmol by atom numbers; sort it later.
         hvec_kmol.resize(n);
         for (int i = 0; i < n; ++i) {
            hvec_kmol[i] = i;
         }

         // sort hvec_kmol
         // put all of the "rattle atoms" in front of the "free atoms"
         auto hvec_kmol_sorter = [&](int ai, int aj) {
            // assign big numbers for atom i and j; if atom i and j are found in
            // constraints, assign smaller numbers to them.
            int ordi = n + ai, ordj = n + aj;
            auto end = hc_dict.end();
            auto it = hc_dict.find(ai);
            auto jt = hc_dict.find(aj);
            if (it != end) {
               ordi = it->second;
            }
            if (jt != end) {
               ordj = jt->second;
            }
            // if atom i and j are in the same "rattle molecule",
            // sort by their atom numbers
            if (ordi == ordj) {
               return ai < aj;
            } else {
               return ordi < ordj;
            }
         };
         std::sort(hvec_kmol.begin(), hvec_kmol.end(), hvec_kmol_sorter);

         hvec_molec.resize(n);
         int current_mol = -1;
         rattle_dmol.totmass = 0.0;
         for (int i = 0; i < n; ++i) {
            int k = hvec_kmol[i]; // the i-th atom in hvec_kmol
            double kmass = atomid::mass[k];
            rattle_dmol.totmass += kmass;
            auto kt = hc_dict.find(k);
            if (kt != hc_dict.end()) {
               // atom k is a "rattle atom"
               int k_mol = kt->second;
               if (k_mol != current_mol) {
                  // we are the first atom of a new "rattle molecule"
                  hvec_molmass.push_back(kmass);
                  current_mol = k_mol;
                  hvec_imol.push_back(i);
                  hvec_imol.push_back(i + 1);
               } else {
                  // we are still at the same "rattle molecule"
                  hvec_molmass.back() += kmass;
                  auto& iend = hvec_imol.back();
                  ++iend;
               }
            } else {
               // atom k is not a "rattle atom", and we will not meet another
               // "rattle atom".
               hvec_molmass.push_back(kmass);
               ++current_mol;
               hvec_imol.push_back(i);
               hvec_imol.push_back(i + 1);
            }
            hvec_molec[k] = current_mol;
         }
         hvec_massfrac.resize(n);
         for (int i = 0; i < n; ++i) {
            double imass = mass[i];
            int im = hvec_molec[i];
            double mmass = hvec_molmass[im];
            hvec_massfrac[i] = imass / mmass;
         }

         // allocate
         if (rc_flag & calc::md) {
            // Actually these arrays are only used for NPT RATTLE.
            darray::allocate(n, &ratcom_x, &ratcom_y, &ratcom_z, &ratcom_px,
                             &ratcom_py, &ratcom_pz, &ratcom_gx, &ratcom_gy,
                             &ratcom_gz, &ratcom_massfrac);
         }

         int nrmol = hvec_molmass.size();
         int nmol2 = hvec_imol.size();
         assert(2 * nrmol == nmol2 && "Error in RATTLE setup.");

         rattle_dmol.nmol = nrmol;
         darray::allocate(nrmol, &rattle_dmol.imol, &rattle_dmol.molmass);
         darray::allocate(n, &rattle_dmol.kmol, &rattle_dmol.molecule);
         darray::copyin(g::q0, nrmol, rattle_dmol.imol, hvec_imol.data());
         darray::copyin(g::q0, nrmol, rattle_dmol.molmass, hvec_molmass.data());
         darray::copyin(g::q0, n, rattle_dmol.kmol, hvec_kmol.data());
         darray::copyin(g::q0, n, rattle_dmol.molecule, hvec_molec.data());
         darray::copyin(g::q0, n, ratcom_massfrac, hvec_massfrac.data());
         wait_for(g::q0);
      }


      rateps = freeze::rateps;


      // find water-like constraints in hc_mols
      // find -CH, -CH2, -CH3
      std::vector<int> veciwater;
      std::vector<pos_prec> veckwater;
      std::vector<int> vecich;
      std::vector<pos_prec> veckch;
      std::vector<int> vecich2;
      std::vector<pos_prec> veckch2;
      std::vector<int> vecich3;
      std::vector<pos_prec> veckch3;
      for (auto& it : hc_mols) {
         int a, b, c, d;
         pos_prec ab, ac, bc, ad;
         if (HCMol_is_water(it, a, b, c, ab, ac, bc)) {
            veciwater.push_back(a);
            veciwater.push_back(b);
            veciwater.push_back(c);
            veckwater.push_back(ab);
            veckwater.push_back(ac);
            veckwater.push_back(bc);
            it.clear();
         } else if (HCMol_is_ch(it, a, b, ab)) {
            vecich.push_back(a);
            vecich.push_back(b);
            veckch.push_back(ab);
            it.clear();
         } else if (HCMol_is_ch2(it, a, b, c, ab, ac)) {
            vecich2.push_back(a);
            vecich2.push_back(b);
            vecich2.push_back(c);
            veckch2.push_back(ab);
            veckch2.push_back(ac);
            if (TINKER_CUDART and pltfm_config & CU_PLTFM)
               it.clear();
         } else if (HCMol_is_ch3(it, a, b, c, d, ab, ac, ad)) {
            vecich3.push_back(a);
            vecich3.push_back(b);
            vecich3.push_back(c);
            vecich3.push_back(d);
            veckch3.push_back(ab);
            veckch3.push_back(ac);
            veckch3.push_back(ad);
            if (TINKER_CUDART and pltfm_config & CU_PLTFM)
               it.clear();
         }
      }
      assert(veciwater.size() % 3 == 0);
      nratwt = veciwater.size() / 3;
      darray::allocate(nratwt, &iratwt, &kratwt);
      darray::copyin(g::q0, nratwt, iratwt, veciwater.data());
      darray::copyin(g::q0, nratwt, kratwt, veckwater.data());
      wait_for(g::q0);
      assert(vecich.size() % 2 == 0);
      nratch = vecich.size() / 2;
      darray::allocate(nratch, &iratch, &kratch);
      darray::copyin(g::q0, nratch, iratch, vecich.data());
      darray::copyin(g::q0, nratch, kratch, veckch.data());
      wait_for(g::q0);
      assert(vecich2.size() % 3 == 0);
      nratch2 = vecich2.size() / 3;
      darray::allocate(nratch2, &iratch2, &kratch2);
      darray::copyin(g::q0, nratch2, iratch2, vecich2.data());
      darray::copyin(g::q0, nratch2, kratch2, veckch2.data());
      wait_for(g::q0);
      assert(vecich3.size() % 4 == 0);
      nratch3 = vecich3.size() / 4;
      darray::allocate(nratch3, &iratch3, &kratch3);
      darray::copyin(g::q0, nratch3, iratch3, vecich3.data());
      darray::copyin(g::q0, nratch3, kratch3, veckch3.data());
      wait_for(g::q0);


      // erase water-like and methyl-like constraints in hc_mols
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
      darray::copyin(g::q0, nrat, irat, iratn.data());
      darray::copyin(g::q0, nrat, krat, kr.data());
      darray::copyin(g::q0, nratmol, iratmol, iratm.data());
      wait_for(g::q0);


      darray::allocate(n, &rattle_xold, &rattle_yold, &rattle_zold);


      hc_mols.clear();
      hc_dict.clear();
   }
}


void rattle(time_prec dt, const pos_prec* xold, const pos_prec* yold,
            const pos_prec* zold)
{
   rattle_settle_acc(dt, xold, yold, zold);
   rattle_ch_acc(dt, xold, yold, zold);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      rattle_methyl_cu(dt, xold, yold, zold);
#endif
   rattle_acc(dt, xold, yold, zold);
}


void rattle2(time_prec dt, bool do_v)
{
   if (do_v) {
      darray::zero(g::q0, buffer_size(), vir_buf);
   }


   rattle2_settle_acc(dt, do_v);
   rattle2_ch_acc(dt, do_v);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      rattle2_methyl_cu(dt, do_v);
#endif
   rattle2_acc(dt, do_v);


   if (do_v) {
      virial_prec v[9];
      virial_reduce(v, vir_buf);
      for (int iv = 0; iv < 9; ++iv) {
         vir[iv] += v[iv];
      }
   }
}


void shake(time_prec dt, pos_prec* xnew, pos_prec* ynew, pos_prec* znew,
           const pos_prec* xold, const pos_prec* yold, const pos_prec* zold)
{
   shake_settle_acc(dt, xnew, ynew, znew, xold, yold, zold);
   shake_ch_acc(dt, xnew, ynew, znew, xold, yold, zold);
#if TINKER_CUDART
   if (pltfm_config & CU_PLTFM)
      shake_methyl_cu(dt, xnew, ynew, znew, xold, yold, zold);
#endif
   shake_acc(dt, xnew, ynew, znew, xold, yold, zold);
}
}
