#include "add.h"
#include "echarge.h"
#include "launch.h"
#include "mdegv.h"
#include "mdpq.h"
#include "named_struct.h"
#include "pme.h"
#include "seq_image.h"
#include "seq_pair_charge.h"
#include "spatial.h"
#include "switch.h"


TINKER_NAMESPACE_BEGIN
#define ECHARGE_ARGS                                                           \
   size_t bufsize, count_buffer restrict nec, energy_buffer restrict ec,       \
      virial_buffer restrict vir_ec, grad_prec *restrict gx,                   \
      grad_prec *restrict gy, grad_prec *restrict gz, TINKER_IMAGE_PARAMS,     \
      real off, real ebuffer, real f, const real *restrict pchg


template <class Ver, class ETYP>
__global__
void echarge_cu1(ECHARGE_ARGS, const Spatial::SortedAtom* restrict sorted,
                 int niak, const int* restrict iak, const int* restrict lst,
                 int n, real aewald)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   const int iwarp = (threadIdx.x + blockIdx.x * blockDim.x) / WARP_SIZE;
   const int nwarp = blockDim.x * gridDim.x / WARP_SIZE;
   const int ilane = threadIdx.x & (WARP_SIZE - 1);
   const int offset = (threadIdx.x + blockIdx.x * blockDim.x) & (bufsize - 1);


   struct Data
   {
      real x, y, z, chg;
      real frcx, frcy, frcz;
      real padding_;
   };
   __shared__ Data data[BLOCK_DIM];


   const real off2 = off * off;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      int ctl;
      if CONSTEXPR (do_a) {
         ctl = 0;
      }
      real etl;
      if CONSTEXPR (do_e) {
         etl = 0;
      }
      real vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz;
      if CONSTEXPR (do_v) {
         vtlxx = 0;
         vtlxy = 0;
         vtlxz = 0;
         vtlyy = 0;
         vtlyz = 0;
         vtlzz = 0;
      }


      Data idat;
      if CONSTEXPR (do_g) {
         idat.frcx = 0;
         idat.frcy = 0;
         idat.frcz = 0;
      }
      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      idat.x = sorted[atomi].x;
      idat.y = sorted[atomi].y;
      idat.z = sorted[atomi].z;
      int i = sorted[atomi].unsorted;
      idat.chg = pchg[i];


      if CONSTEXPR (do_g) {
         data[threadIdx.x].frcx = 0;
         data[threadIdx.x].frcy = 0;
         data[threadIdx.x].frcz = 0;
      }
      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].x = sorted[shatomk].x;
      data[threadIdx.x].y = sorted[shatomk].y;
      data[threadIdx.x].z = sorted[shatomk].z;
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].chg = pchg[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real xr = idat.x - data[klane].x;
         real yr = idat.y - data[klane].y;
         real zr = idat.z - data[klane].z;


         real r2 = image2(xr, yr, zr);
         if (atomi < atomk && r2 <= off2) {
            real r = REAL_SQRT(r2);


            MAYBE_UNUSED real grdx, grdy, grdz;
            if CONSTEXPR (do_g)
               grdx = grdy = grdz = 0;


            if CONSTEXPR (eq<ETYP, EWALD>()) {
               pair_charge<Ver, EWALD>(r, xr, yr, zr, 1, idat.chg,
                                       data[klane].chg, ebuffer, f, aewald, //
                                       grdx, grdy, grdz, ctl, etl, vtlxx, vtlxy,
                                       vtlxz, vtlyy, vtlyz, vtlzz);
            }


            if CONSTEXPR (do_g) {
               idat.frcx += grdx;
               idat.frcy += grdy;
               idat.frcz += grdz;
               data[klane].frcx -= grdx;
               data[klane].frcy -= grdy;
               data[klane].frcz -= grdz;
            }
         } // end if (include)
      }


      if CONSTEXPR (do_a)
         atomic_add(ctl, nec, offset);
      if CONSTEXPR (do_e)
         atomic_add(etl, ec, offset);
      if CONSTEXPR (do_g) {
         atomic_add(idat.frcx, &gx[i]);
         atomic_add(idat.frcy, &gy[i]);
         atomic_add(idat.frcz, &gz[i]);
         atomic_add(data[threadIdx.x].frcx, &gx[shk]);
         atomic_add(data[threadIdx.x].frcy, &gy[shk]);
         atomic_add(data[threadIdx.x].frcz, &gz[shk]);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz, vir_ec, offset);
   } // end for (iw)
}


template <class Ver, class ETYP>
__global__
void echarge_cu2(ECHARGE_ARGS, const real* restrict x, const real* restrict y,
                 const real* restrict z, int ncexclude,
                 const int (*restrict cexclude)[2],
                 const real* restrict cexclude_scale)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < ncexclude;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);


      int i = cexclude[ii][0];
      int k = cexclude[ii][1];
      real cscale = cexclude_scale[ii];


      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real ci = pchg[i];


      real xr = xi - x[k];
      real yr = yi - y[k];
      real zr = zi - z[k];
      real ck = pchg[k];


      real r2 = image2(xr, yr, zr);
      real off2 = off * off;
      if (r2 <= off2) {
         MAYBE_UNUSED int ctl;
         MAYBE_UNUSED real e, grdx, grdy, grdz, vxx, vxy, vxz, vyy, vyz, vzz;
         if CONSTEXPR (do_a) {
            ctl = 0;
         }
         if CONSTEXPR (do_e) {
            e = 0;
         }
         if CONSTEXPR (do_g) {
            grdx = 0;
            grdy = 0;
            grdz = 0;
         }
         if CONSTEXPR (do_v) {
            vxx = 0;
            vxy = 0;
            vxz = 0;
            vyy = 0;
            vyz = 0;
            vzz = 0;
         }


         real r = REAL_SQRT(r2);
         pair_charge<Ver, ETYP>(r, xr, yr, zr, cscale, ci, ck, ebuffer, f, 0, //
                                grdx, grdy, grdz, ctl, e, vxx, vxy, vxz, vyy,
                                vyz, vzz);
         if (e != 0) {
            if CONSTEXPR (do_a) {
               atomic_add(ctl, nec, offset);
            }
            if CONSTEXPR (do_e) {
               atomic_add(e, ec, offset);
            }
            if CONSTEXPR (do_g) {
               atomic_add(grdx, gx, i);
               atomic_add(grdy, gy, i);
               atomic_add(grdz, gz, i);
               atomic_add(-grdx, gx, k);
               atomic_add(-grdy, gy, k);
               atomic_add(-grdz, gz, k);
            }
            if CONSTEXPR (do_v) {
               atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_ec, offset);
            }
         }
      } // end if (include)
   }
}


template <class Ver, class ETYP>
void echarge_cu()
{
   const auto& st = *cspatial_unit;
   real off;
   if CONSTEXPR (eq<ETYP, EWALD>())
      off = switch_off(switch_ewald);
   else
      off = switch_off(switch_charge);
   auto bufsize = buffer_size();


   const real f = electric / dielec;
   real aewald = 0;
   if CONSTEXPR (eq<ETYP, EWALD>()) {
      PMEUnit pu = epme_unit;
      aewald = pu->aewald;
   }


   if CONSTEXPR (eq<ETYP, EWALD>()) {
      if (st.niak > 0) {
         auto ker1 = echarge_cu1<Ver, EWALD>;
         launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                    bufsize, nec, ec, vir_ec, gx, gy, gz, TINKER_IMAGE_ARGS,
                    off, ebuffer, f, pchg, //
                    st.sorted, st.niak, st.iak, st.lst, n, aewald);
      }
      if (ncexclude > 0) {
         auto ker2 = echarge_cu2<Ver, NON_EWALD>;
         launch_k1s(nonblk, ncexclude, ker2, //
                    bufsize, nec, ec, vir_ec, gx, gy, gz, TINKER_IMAGE_ARGS,
                    off, ebuffer, f, pchg, //
                    x, y, z, ncexclude, cexclude, cexclude_scale);
      }
   } else if CONSTEXPR (eq<ETYP, NON_EWALD_TAPER>()) {
   }
}


void echarge_ewald_real_cu(int vers)
{
   if (vers == calc::v0)
      echarge_cu<calc::V0, EWALD>();
   else if (vers == calc::v1)
      echarge_cu<calc::V1, EWALD>();
   else if (vers == calc::v3)
      echarge_cu<calc::V3, EWALD>();
   else if (vers == calc::v4)
      echarge_cu<calc::V4, EWALD>();
   else if (vers == calc::v5)
      echarge_cu<calc::V5, EWALD>();
   else if (vers == calc::v6)
      echarge_cu<calc::V6, EWALD>();
}


//====================================================================//


template <class Ver>
__global__
void echarge_cu3(size_t bufsize, count_buffer restrict nec,
                 energy_buffer restrict ec, const real* restrict pchg, real f,
                 real aewald, int n)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;


   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < n;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);
      real chgi = pchg[ii];
      if (chgi == 0)
         continue;


      // self energy, tinfoil
      if CONSTEXPR (do_e) {
         real fs = -f * aewald * REAL_RECIP(sqrtpi);
         real e = fs * chgi * chgi;
         atomic_add(e, ec, offset);
         if (do_a) {
            atomic_add(1, nec, offset);
         }
      }


      if CONSTEXPR (do_g) {
      }
   }
}


template <class Ver>
void echarge_fphi_self_cu()
{
   auto bufsize = buffer_size();
   real f = electric / dielec;
   real aewald = epme_unit->aewald;


   auto ker = echarge_cu3<Ver>;
   launch_k2s(nonblk, PME_BLOCKDIM, n,
              ker, //
              bufsize, nec, ec, pchg, f, aewald, n);
}


void echarge_ewald_fphi_self_cu(int vers)
{
   if (vers == calc::v0)
      echarge_fphi_self_cu<calc::V0>();
   else if (vers == calc::v1)
      echarge_fphi_self_cu<calc::V1>();
   else if (vers == calc::v3)
      echarge_fphi_self_cu<calc::V3>();
   else if (vers == calc::v4)
      echarge_fphi_self_cu<calc::V4>();
   else if (vers == calc::v5)
      echarge_fphi_self_cu<calc::V5>();
   else if (vers == calc::v6)
      echarge_fphi_self_cu<calc::V6>();
}
TINKER_NAMESPACE_END
