#include "add.h"
#include "evdw.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "named_struct.h"
#include "seq_pair_lj.h"
#include "seq_switch.h"
#include "spatial.h"
#include "switch.h"


TINKER_NAMESPACE_BEGIN
#define LJ_ARGS                                                                \
   size_t bufsize, count_buffer restrict nev, energy_buffer restrict ev,       \
      virial_buffer restrict vir_ev, grad_prec *restrict gx,                   \
      grad_prec *restrict gy, grad_prec *restrict gz, TINKER_IMAGE_PARAMS,     \
      int njvdw, const int *restrict jvdw, const real *restrict radmin,        \
      const real *restrict epsilon, real cut, real off


template <class Ver>
__global__
void elj_cu1(LJ_ARGS, int n, const Spatial::SortedAtom* restrict sorted,
             int niak, const int* restrict iak, const int* restrict lst)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int iwarp = ithread / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = ithread & (bufsize - 1);


   // thread local variables
   using ebuf_prec = energy_buffer_traits::type;
   using vbuf_prec = virial_buffer_traits::type;
   MAYBE_UNUSED int ctl;
   MAYBE_UNUSED ebuf_prec etl;
   MAYBE_UNUSED grad_prec gxi, gyi, gzi, gxk, gyk, gzk;
   MAYBE_UNUSED vbuf_prec vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      if CONSTEXPR (do_a)
         ctl = 0;
      if CONSTEXPR (do_e)
         etl = 0;
      if CONSTEXPR (do_g) {
         gxi = 0;
         gyi = 0;
         gzi = 0;
         gxk = 0;
         gyk = 0;
         gzk = 0;
      }
      if CONSTEXPR (do_v) {
         vtlxx = 0;
         vtlyx = 0;
         vtlzx = 0;
         vtlyy = 0;
         vtlzy = 0;
         vtlzz = 0;
      }


      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      real xi = sorted[atomi].x;
      real yi = sorted[atomi].y;
      real zi = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      int it = jvdw[i];


      int shatomk = lst[iw * WARP_SIZE + ilane];
      real shx = sorted[shatomk].x;
      real shy = sorted[shatomk].y;
      real shz = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      int shkt = jvdw[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = xi - __shfl_sync(ALL_LANES, shx, srclane);
         real yr = yi - __shfl_sync(ALL_LANES, shy, srclane);
         real zr = zi - __shfl_sync(ALL_LANES, shz, srclane);
         int kt = __shfl_sync(ALL_LANES, shkt, srclane);


         MAYBE_UNUSED real dedx = 0, dedy = 0, dedz = 0;


         real rik2 = image2(xr, yr, zr);
         if (atomi < atomk && rik2 <= off2) {
            real rik = REAL_SQRT(rik2);
            real rv = radmin[it * njvdw + kt];
            real eps = epsilon[it * njvdw + kt];


            MAYBE_UNUSED real e, de;
            pair_lj<do_g>(rik, rik2, rv, eps, 1, e, de);


            if (rik2 > cut2) {
               real taper, dtaper;
               switch_taper5<do_g>(rik, cut, off, taper, dtaper);
               if CONSTEXPR (do_g)
                  de = e * dtaper + de * taper;
               if CONSTEXPR (do_e)
                  e = e * taper;
            }


            if CONSTEXPR (do_a)
               if (e != 0)
                  ctl += 1;
            if CONSTEXPR (do_e)
               etl += to_cu<ebuf_prec>(e);
            if CONSTEXPR (do_g) {
               de *= REAL_RECIP(rik);
               dedx = de * xr;
               dedy = de * yr;
               dedz = de * zr;
               if CONSTEXPR (do_v) {
                  vtlxx += to_cu<vbuf_prec>(xr * dedx);
                  vtlyx += to_cu<vbuf_prec>(yr * dedx);
                  vtlzx += to_cu<vbuf_prec>(zr * dedx);
                  vtlyy += to_cu<vbuf_prec>(yr * dedy);
                  vtlzy += to_cu<vbuf_prec>(zr * dedy);
                  vtlzz += to_cu<vbuf_prec>(zr * dedz);
               }
            }
         } // end if (include)


         if CONSTEXPR (do_g) {
            int dstlane = (ilane + WARP_SIZE - j) & (WARP_SIZE - 1);
            gxi += to_cu<grad_prec>(dedx);
            gyi += to_cu<grad_prec>(dedy);
            gzi += to_cu<grad_prec>(dedz);
            gxk -= to_cu<grad_prec>(__shfl_sync(ALL_LANES, dedx, dstlane));
            gyk -= to_cu<grad_prec>(__shfl_sync(ALL_LANES, dedy, dstlane));
            gzk -= to_cu<grad_prec>(__shfl_sync(ALL_LANES, dedz, dstlane));
         }
         __syncwarp();
      }


      if CONSTEXPR (do_a)
         atomic_add(ctl, nev, offset);
      if CONSTEXPR (do_e)
         atomic_add(etl, ev, offset);
      if CONSTEXPR (do_g) {
         atomic_add(gxi, gx, i);
         atomic_add(gyi, gy, i);
         atomic_add(gzi, gz, i);
         atomic_add(gxk, gx, shk);
         atomic_add(gyk, gy, shk);
         atomic_add(gzk, gz, shk);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vir_ev, offset);
   } // end for (iw)
}


template <class Ver>
__global__
void elj_cu2(LJ_ARGS, const real* restrict x, const real* restrict y,
             const real* restrict z, int nvexclude,
             const int (*restrict vexclude)[2],
             const real* restrict vexclude_scale)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nvexclude;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = vexclude[ii][0];
      int k = vexclude[ii][1];
      real vscale = vexclude_scale[ii];


      int it = jvdw[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      int kt = jvdw[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];


      real rik2 = image2(xr, yr, zr);
      if (rik2 <= off2) {
         real rik = REAL_SQRT(rik2);
         real rv = radmin[it * njvdw + kt];
         real eps = epsilon[it * njvdw + kt];


         MAYBE_UNUSED real e, de;
         pair_lj<do_g>(rik, rik2, rv, eps, vscale, e, de);
         if (rik2 > cut2) {
            real taper, dtaper;
            switch_taper5<do_g>(rik, cut, off, taper, dtaper);
            if CONSTEXPR (do_g)
               de = e * dtaper + de * taper;
            if CONSTEXPR (do_e)
               e = e * taper;
         }


         if CONSTEXPR (do_a)
            if (vscale == -1 && e != 0)
               atomic_add(-1, nev, offset);
         if CONSTEXPR (do_e)
            atomic_add(e, ev, offset);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(rik);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;
            atomic_add(dedx, gx, i);
            atomic_add(dedy, gy, i);
            atomic_add(dedz, gz, i);
            atomic_add(-dedx, gx, k);
            atomic_add(-dedy, gy, k);
            atomic_add(-dedz, gz, k);
            if CONSTEXPR (do_v) {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;
               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ev, offset);
            }
         }
      } // end if (include)
      __syncwarp();
   }
}


template <class Ver>
__global__
void elj_cu3(LJ_ARGS, const real* restrict x, const real* restrict y,
             const real* restrict z, real v4scale, int nvdw14,
             const int (*restrict vdw14ik)[2], const real* restrict radmin4,
             const real* restrict epsilon4)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const real cut2 = cut * cut;
   const real off2 = off * off;
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nvdw14;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = vdw14ik[ii][0];
      int k = vdw14ik[ii][1];


      int it = jvdw[i];
      real xi = x[i];
      real yi = y[i];
      real zi = z[i];


      int kt = jvdw[k];
      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];


      real rik2 = image2(xr, yr, zr);
      if (rik2 <= off2) {
         real rik = REAL_SQRT(rik2);
         real rv = radmin[it * njvdw + kt];
         real eps = epsilon[it * njvdw + kt];
         real rv4 = radmin4[it * njvdw + kt];
         real eps4 = epsilon4[it * njvdw + kt];


         MAYBE_UNUSED real e, de, e4, de4;
         pair_lj<do_g>(rik, rik2, rv, eps, v4scale, e, de);
         pair_lj<do_g>(rik, rik2, rv4, eps4, v4scale, e4, de4);
         e = e4 - e;
         if CONSTEXPR (do_g)
            de = de4 - de;


         if (rik2 > cut2) {
            real taper, dtaper;
            switch_taper5<do_g>(rik, cut, off, taper, dtaper);
            if CONSTEXPR (do_g)
               de = e * dtaper + de * taper;
            if CONSTEXPR (do_e)
               e = e * taper;
         }


         // if CONSTEXPR (do_a) {}
         if CONSTEXPR (do_e)
            atomic_add(e, ev, offset);
         if CONSTEXPR (do_g) {
            de *= REAL_RECIP(rik);
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;
            atomic_add(dedx, gx, i);
            atomic_add(dedy, gy, i);
            atomic_add(dedz, gz, i);
            atomic_add(-dedx, gx, k);
            atomic_add(-dedy, gy, k);
            atomic_add(-dedz, gz, k);
            if CONSTEXPR (do_v) {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;
               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_ev, offset);
            }
         }
      } // end if (include)
      __syncwarp();
   }
}


template <class Ver>
void elj_cu4()
{
   const auto& st = *cspatial_unit;
   const real cut = switch_cut(switch_vdw);
   const real off = switch_off(switch_vdw);


   auto bufsize = buffer_size();


   if (st.niak > 0)
      launch_k1s(nonblk, WARP_SIZE * st.niak, elj_cu1<Ver>, bufsize, nev, ev,
                 vir_ev, gx, gy, gz, TINKER_IMAGE_ARGS, njvdw, jvdw, radmin,
                 epsilon, cut, off, n, st.sorted, st.niak, st.iak, st.lst);


   if (nvexclude > 0)
      launch_k1s(nonblk, nvexclude, elj_cu2<Ver>, bufsize, nev, ev, vir_ev, gx,
                 gy, gz, TINKER_IMAGE_ARGS, njvdw, jvdw, radmin, epsilon, cut,
                 off, x, y, z, nvexclude, vexclude, vexclude_scale);


   if (nvdw14 > 0)
      launch_k1s(nonblk, nvdw14, elj_cu3<Ver>, bufsize, nev, ev, vir_ev, gx, gy,
                 gz, TINKER_IMAGE_ARGS, njvdw, jvdw, radmin, epsilon, cut, off,
                 x, y, z, v4scale, nvdw14, vdw14ik, radmin4, epsilon4);
}


void elj_cu(int vers)
{
   if (vers == calc::v0)
      elj_cu4<calc::V0>();
   else if (vers == calc::v1)
      elj_cu4<calc::V1>();
   else if (vers == calc::v3)
      elj_cu4<calc::V3>();
   else if (vers == calc::v4)
      elj_cu4<calc::V4>();
   else if (vers == calc::v5)
      elj_cu4<calc::V5>();
   else if (vers == calc::v6)
      elj_cu4<calc::V6>();
}
TINKER_NAMESPACE_END
