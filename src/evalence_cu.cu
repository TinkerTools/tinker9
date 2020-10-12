#include "add.h"
#include "energy.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "potent.h"
#include "seq_torsion.h"
#include "tool/gpu_card.h"
#include "tool/host_zero.h"


namespace tinker {
template <class Ver>
__global__
void evalence_cu1(
   // torsion
   energy_buffer restrict et, virial_buffer restrict vir_et,
   grad_prec* restrict detx, grad_prec* restrict dety, grad_prec* restrict detz,
   real torsunit, int ntors, const int (*restrict itors)[4],
   const real (*restrict tors1)[4], const real (*restrict tors2)[4],
   const real (*restrict tors3)[4], const real (*restrict tors4)[4],
   const real (*restrict tors5)[4], const real (*restrict tors6)[4],
   // xyz
   const real* restrict x, const real* restrict y, const real* restrict z)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int stride = blockDim.x * gridDim.x;


   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec e0t; // torsion
   if CONSTEXPR (do_e) {
      e0t = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec v0txx, v0tyx, v0tzx, v0tyy, v0tzy, v0tzz;
   if CONSTEXPR (do_v) {
      v0txx = 0, v0tyx = 0, v0tzx = 0, v0tyy = 0, v0tzy = 0, v0tzz = 0;
   }


   // torsion
   for (int i = ithread; i < ntors; i += stride) {
      int ia, ib, ic, id;
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      real dedxia, dedyia, dedzia;
      real dedxib, dedyib, dedzib;
      real dedxic, dedyic, dedzic;
      real dedxid, dedyid, dedzid;
      dk_tors<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

                   dedxia, dedyia, dedzia, dedxib, dedyib, dedzib, dedxic,
                   dedyic, dedzic, dedxid, dedyid, dedzid,

                   torsunit, i, ia, ib, ic, id, itors,

                   tors1, tors2, tors3, tors4, tors5, tors6, x, y, z);
      if CONSTEXPR (do_e) {
         e0t += cvt_to<ebuf_prec>(e);
      }
      if CONSTEXPR (do_g) {
         atomic_add(dedxia, detx, ia);
         atomic_add(dedyia, dety, ia);
         atomic_add(dedzia, detz, ia);
         atomic_add(dedxib, detx, ib);
         atomic_add(dedyib, dety, ib);
         atomic_add(dedzib, detz, ib);
         atomic_add(dedxic, detx, ic);
         atomic_add(dedyic, dety, ic);
         atomic_add(dedzic, detz, ic);
         atomic_add(dedxid, detx, id);
         atomic_add(dedyid, dety, id);
         atomic_add(dedzid, detz, id);
      }
      if CONSTEXPR (do_v) {
         v0txx += cvt_to<vbuf_prec>(vxx);
         v0tyx += cvt_to<vbuf_prec>(vyx);
         v0tzx += cvt_to<vbuf_prec>(vzx);
         v0tyy += cvt_to<vbuf_prec>(vyy);
         v0tzy += cvt_to<vbuf_prec>(vzy);
         v0tzz += cvt_to<vbuf_prec>(vzz);
      }
   }
   if CONSTEXPR (do_e and do_a) {
      atomic_add(e0t, et, ithread);
   }
   if CONSTEXPR (do_v and do_a) {
      atomic_add(v0txx, v0tyx, v0tzx, v0tyy, v0tzy, v0tzz, vir_et, ithread);
   }


   // total energy and virial
   if CONSTEXPR (do_e and not do_a) {
      ebuf_prec etl = 0;
      etl += e0t;
      atomic_add(etl, et, ithread);
   }
   if CONSTEXPR (do_v and not do_a) {
      vbuf_prec vtlxx = 0, vtlyx = 0, vtlzx = 0;
      vbuf_prec vtlyy = 0, vtlzy = 0, vtlzz = 0;
      vtlxx += v0txx, vtlyx += v0tyx, vtlzx += v0tzx;
      vtlyy += v0tyy, vtlzy += v0tzy, vtlzz += v0tzz;
      atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vir_et, ithread);
   }
}


void evalence_cu2(int vers)
{
#define EVALENCE_ARGS                                                          \
   et, vir_et, detx, dety, detz, torsunit, ntors, itors, tors1, tors2, tors3,  \
      tors4, tors5, tors6, x, y, z


   int ngrid = get_grid_size(BLOCK_DIM);
   if (vers == calc::v0)
      evalence_cu1<calc::V0><<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
   else if (vers == calc::v1)
      evalence_cu1<calc::V1><<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
   else if (vers == calc::v3)
      evalence_cu1<calc::V3><<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
   else if (vers == calc::v4)
      evalence_cu1<calc::V4><<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
   else if (vers == calc::v5)
      evalence_cu1<calc::V5><<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
   else if (vers == calc::v6)
      evalence_cu1<calc::V6><<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
}


void evalence_cu(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   bool flag_torsion = use_potent(torsion_term);


   size_t bsize = buffer_size();
   if (rc_a and flag_torsion) {
      host_zero(energy_et, virial_et);
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, et);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_et);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, detx, dety, detz);
   }


   if (flag_torsion) {
      evalence_cu2(vers);
   }


   if (rc_a and flag_torsion) {
      if (do_e) {
         energy_et = energy_reduce(et);
         energy_valence += energy_et;
      }
      if (do_v) {
         virial_reduce(virial_et, vir_et);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_et[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, detx, dety, detz);
   }
}
}
