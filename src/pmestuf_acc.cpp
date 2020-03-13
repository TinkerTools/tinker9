#include "add.h"
#include "box.h"
#include "elec.h"
#include "mdpq.h"
#include "named_struct.h"
#include "pmestuf.h"
#include "seq_bsplgen.h"


TINKER_NAMESPACE_BEGIN
template <class T>
void grid_put_acc(PMEUnit pme_u, real* ptr1, real* ptr2)
{
   auto& st = *pme_u;
   auto* qgrid = st.qgrid;


   MAYBE_UNUSED const real* pchg = ptr1;
   MAYBE_UNUSED const real(*fmp)[10] = (real(*)[10])ptr1;
   MAYBE_UNUSED const real(*fuind)[3] = (real(*)[3])ptr1;
   MAYBE_UNUSED const real(*fuinp)[3] = (real(*)[3])ptr2;


   const int nfft1 = st.nfft1;
   const int nfft2 = st.nfft2;
   const int nfft3 = st.nfft3;
   const int bsorder = st.bsorder;
   assert(bsorder <= 5);


   darray::zero(PROCEED_NEW_Q, 2 * nfft1 * nfft2 * nfft3, st.qgrid);


   #pragma acc parallel loop independent async\
               deviceptr(pchg,fmp,fuind,fuinp,x,y,z,qgrid)
   for (int i = 0; i < n; ++i) {
      real thetai1[4 * 5];
      real thetai2[4 * 5];
      real thetai3[4 * 5];


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      // map fractional coordinate w from [-0.5 + k, 0.5 + k) to [0, 1)
      // w -> (w + 0.5) - FLOOR(w + 0.5)
      // see also subroutine bspline_fill in pmestuf.f


      real w1 = xi * recipa.x + yi * recipa.y + zi * recipa.z;
      w1 = w1 + 0.5f - REAL_FLOOR(w1 + 0.5f);
      real fr1 = nfft1 * w1;
      int igrid1 = REAL_FLOOR(fr1);
      w1 = fr1 - igrid1;


      real w2 = xi * recipb.x + yi * recipb.y + zi * recipb.z;
      w2 = w2 + 0.5f - REAL_FLOOR(w2 + 0.5f);
      real fr2 = nfft2 * w2;
      int igrid2 = REAL_FLOOR(fr2);
      w2 = fr2 - igrid2;


      real w3 = xi * recipc.x + yi * recipc.y + zi * recipc.z;
      w3 = w3 + 0.5f - REAL_FLOOR(w3 + 0.5f);
      real fr3 = nfft3 * w3;
      int igrid3 = REAL_FLOOR(fr3);
      w3 = fr3 - igrid3;


      igrid1 = igrid1 - bsorder + 1;
      igrid2 = igrid2 - bsorder + 1;
      igrid3 = igrid3 - bsorder + 1;
      igrid1 += (igrid1 < 0 ? nfft1 : 0);
      igrid2 += (igrid2 < 0 ? nfft2 : 0);
      igrid3 += (igrid3 < 0 ? nfft3 : 0);


      if CONSTEXPR (eq<T, PCHG>()) {
         real pchgi = pchg[i];
         if (pchgi == 0)
            continue;


         bsplgen<1>(w1, thetai1, bsorder);
         bsplgen<1>(w2, thetai2, bsorder);
         bsplgen<1>(w3, thetai3, bsorder);


         #pragma acc loop seq
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real v0 = thetai3[4 * iz] * pchgi;
            #pragma acc loop seq
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real u0 = thetai2[4 * iy] * v0;
               #pragma acc loop seq
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  int index = xbase + ybase + zbase;
                  real term = thetai1[4 * ix] * u0;
                  atomic_add(term, qgrid, 2 * index);
               }
            } // end for (int iy)
         }
      } // end if (grid_pchg)


      if CONSTEXPR (eq<T, MPOLE>()) {
         bsplgen<3>(w1, thetai1, bsorder);
         bsplgen<3>(w2, thetai2, bsorder);
         bsplgen<3>(w3, thetai3, bsorder);


         real fmpi0 = fmp[i][mpl_pme_0];
         real fmpix = fmp[i][mpl_pme_x];
         real fmpiy = fmp[i][mpl_pme_y];
         real fmpiz = fmp[i][mpl_pme_z];
         real fmpixx = fmp[i][mpl_pme_xx];
         real fmpiyy = fmp[i][mpl_pme_yy];
         real fmpizz = fmp[i][mpl_pme_zz];
         real fmpixy = fmp[i][mpl_pme_xy];
         real fmpixz = fmp[i][mpl_pme_xz];
         real fmpiyz = fmp[i][mpl_pme_yz];
         #pragma acc loop seq
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real v0 = thetai3[4 * iz];
            real v1 = thetai3[4 * iz + 1];
            real v2 = thetai3[4 * iz + 2];
            #pragma acc loop seq
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real u0 = thetai2[4 * iy];
               real u1 = thetai2[4 * iy + 1];
               real u2 = thetai2[4 * iy + 2];
               // fmp: 0, x, y, z, xx, yy, zz, xy, xz, yz
               //      1, 2, 3, 4,  5,  6,  7,  8,  9, 10
               real term0 = fmpi0 * u0 * v0 + fmpiy * u1 * v0 +
                  fmpiz * u0 * v1 + fmpiyy * u2 * v0 + fmpizz * u0 * v2 +
                  fmpiyz * u1 * v1;
               real term1 =
                  fmpix * u0 * v0 + fmpixy * u1 * v0 + fmpixz * u0 * v1;
               real term2 = fmpixx * u0 * v0;
               #pragma acc loop seq
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  int index = xbase + ybase + zbase;
                  real t0 = thetai1[4 * ix];
                  real t1 = thetai1[4 * ix + 1];
                  real t2 = thetai1[4 * ix + 2];
                  atomic_add(term0 * t0 + term1 * t1 + term2 * t2, qgrid,
                             2 * index);
               }
            } // end for (int iy)
         }
      } // end if (grid_mpole)


      if CONSTEXPR (eq<T, UIND>()) {
         bsplgen<2>(w1, thetai1, bsorder);
         bsplgen<2>(w2, thetai2, bsorder);
         bsplgen<2>(w3, thetai3, bsorder);


         real fuindi0 = fuind[i][0];
         real fuindi1 = fuind[i][1];
         real fuindi2 = fuind[i][2];
         real fuinpi0 = fuinp[i][0];
         real fuinpi1 = fuinp[i][1];
         real fuinpi2 = fuinp[i][2];
         #pragma acc loop seq
         for (int iz = 0; iz < bsorder; ++iz) {
            int zbase = igrid3 + iz;
            zbase -= (zbase >= nfft3 ? nfft3 : 0);
            zbase *= (nfft1 * nfft2);
            real v0 = thetai3[4 * iz];
            real v1 = thetai3[4 * iz + 1];
            #pragma acc loop seq
            for (int iy = 0; iy < bsorder; ++iy) {
               int ybase = igrid2 + iy;
               ybase -= (ybase >= nfft2 ? nfft2 : 0);
               ybase *= nfft1;
               real u0 = thetai2[4 * iy];
               real u1 = thetai2[4 * iy + 1];
               real term01 = fuindi1 * u1 * v0 + fuindi2 * u0 * v1;
               real term11 = fuindi0 * u0 * v0;
               real term02 = fuinpi1 * u1 * v0 + fuinpi2 * u0 * v1;
               real term12 = fuinpi0 * u0 * v0;
               #pragma acc loop seq
               for (int ix = 0; ix < bsorder; ++ix) {
                  int xbase = igrid1 + ix;
                  xbase -= (xbase >= nfft1 ? nfft1 : 0);
                  int index = xbase + ybase + zbase;
                  real t0 = thetai1[4 * ix];
                  real t1 = thetai1[4 * ix + 1];
                  atomic_add(term01 * t0 + term11 * t1, qgrid, 2 * index);
                  atomic_add(term02 * t0 + term12 * t1, qgrid, 2 * index + 1);
               }
            } // end for (int iy)
         }
      } // end if (grid_uind)
   }
}


void grid_pchg_acc(PMEUnit pme_u, real* pchg)
{
   grid_put_acc<PCHG>(pme_u, pchg, nullptr);
}


void grid_mpole_acc(PMEUnit pme_u, real (*fmp)[10])
{
   grid_put_acc<MPOLE>(pme_u, (real*)fmp, nullptr);
}


void grid_uind_acc(PMEUnit pme_u, real (*fuind)[3], real (*fuinp)[3])
{
   grid_put_acc<UIND>(pme_u, (real*)fuind, (real*)fuinp);
}


//====================================================================//


template <bool DO_E, bool DO_V>
void pme_conv_acc1(PMEUnit pme_u, energy_buffer gpu_e, virial_buffer gpu_vir)
{
   auto& st = *pme_u;
   real(*restrict qgrid)[2] = reinterpret_cast<real(*)[2]>(st.qgrid);
   const real* bsmod1 = st.bsmod1;
   const real* bsmod2 = st.bsmod2;
   const real* bsmod3 = st.bsmod3;


   const int nfft1 = st.nfft1;
   const int nfft2 = st.nfft2;
   const int nfft3 = st.nfft3;
   const int nff = nfft1 * nfft2;
   const int ntot = nfft1 * nfft2 * nfft3;


   const real f = electric / dielec;
   real aewald = st.aewald;
   real pterm = pi / aewald;
   pterm *= pterm;
   real box_volume = volbox();


   auto bufsize = buffer_size();
   #pragma acc parallel loop independent async\
               deviceptr(gpu_e,gpu_vir,qgrid,bsmod1,bsmod2,bsmod3)
   for (int i = 0; i < ntot; ++i) {
      if (i == 0) {
         qgrid[0][0] = 0;
         qgrid[0][1] = 0;
         continue;
      }


      int k3 = i / nff;
      int j = i - k3 * nff;
      int k2 = j / nfft1;
      int k1 = j - k2 * nfft1;


      int r1 = (k1 < (nfft1 + 1) / 2) ? k1 : (k1 - nfft1);
      int r2 = (k2 < (nfft2 + 1) / 2) ? k2 : (k2 - nfft2);
      int r3 = (k3 < (nfft3 + 1) / 2) ? k3 : (k3 - nfft3);


      real h1 = recipa.x * r1 + recipb.x * r2 + recipc.x * r3;
      real h2 = recipa.y * r1 + recipb.y * r2 + recipc.y * r3;
      real h3 = recipa.z * r1 + recipb.z * r2 + recipc.z * r3;
      real hsq = h1 * h1 + h2 * h2 + h3 * h3;


      real gridx = qgrid[i][0];
      real gridy = qgrid[i][1];
      real term = -pterm * hsq;
      real expterm = 0;
      if (term > -50) {
         // TODO: if .not. use_bounds; if octahedron; 2/hsq
         real denom =
            hsq * pi * box_volume * bsmod1[k1] * bsmod2[k2] * bsmod3[k3];
         expterm = REAL_EXP(term) / denom;


         if CONSTEXPR (DO_E || DO_V) {
            real struc2 = gridx * gridx + gridy * gridy;
            real eterm = 0.5f * f * expterm * struc2;
            if (DO_E) {
               atomic_add(eterm, gpu_e, i & (bufsize - 1));
            }
            if (DO_V) {
               real vterm = (2 / hsq) * (1 - term) * eterm;

               real vxx = (h1 * h1 * vterm - eterm);
               real vxy = h1 * h2 * vterm;
               real vxz = h1 * h3 * vterm;
               real vyy = (h2 * h2 * vterm - eterm);
               real vyz = h2 * h3 * vterm;
               real vzz = (h3 * h3 * vterm - eterm);
               atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, gpu_vir,
                          i & (bufsize - 1));
            }
         } // end if (e or v)
      }


      // complete the transformation of the PME grid
      qgrid[i][0] = gridx * expterm;
      qgrid[i][1] = gridy * expterm;
   }
}


void pme_conv_acc(PMEUnit pme_u, energy_buffer gpu_e, virial_buffer gpu_vir)
{
   if (gpu_vir == nullptr) {
      if (gpu_e == nullptr) {
         pme_conv_acc1<false, false>(pme_u, nullptr, nullptr);
      } else {
         pme_conv_acc1<true, false>(pme_u, gpu_e, nullptr);
      }
   } else {
      if (gpu_e == nullptr) {
         pme_conv_acc1<false, true>(pme_u, nullptr, gpu_vir);
      } else {
         pme_conv_acc1<true, true>(pme_u, gpu_e, gpu_vir);
      }
   }
}
TINKER_NAMESPACE_END
