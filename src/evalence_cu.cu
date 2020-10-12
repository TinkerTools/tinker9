#include "add.h"
#include "energy.h"
#include "glob.group.h"
#include "glob.molecule.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "potent.h"
#include "seq_geom.h"
#include "seq_torsion.h"
#include "tool/gpu_card.h"
#include "tool/host_zero.h"


namespace tinker {
template <class Ver, bool rc_a>
__global__
void evalence_cu1(
   // etors
   energy_buffer restrict et, virial_buffer restrict vir_et,
   grad_prec* restrict detx, grad_prec* restrict dety, grad_prec* restrict detz,

   real torsunit, int ntors, const int (*restrict itors)[4],
   const real (*restrict tors1)[4], const real (*restrict tors2)[4],
   const real (*restrict tors3)[4], const real (*restrict tors4)[4],
   const real (*restrict tors5)[4], const real (*restrict tors6)[4],

   // geom
   energy_buffer restrict eg, virial_buffer restrict vir_eg,
   grad_prec* restrict degx, grad_prec* restrict degy, grad_prec* restrict degz,

   int ngfix, const int (*restrict igfix)[2], const real (*restrict gfix)[3],

   // total
   energy_buffer restrict ebuf, virial_buffer restrict vbuf,

   // other
   const real* restrict x, const real* restrict y, const real* restrict z,
   const mass_prec* restrict mass, const int* restrict molec,
   const int (*restrict igrp)[2], const int* restrict kgrp,
   const mass_prec* restrict grpmass, TINKER_IMAGE_PARAMS)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int stride = blockDim.x * gridDim.x;


   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec e0t; // etors
   ebuf_prec e0g; // egeom
   if CONSTEXPR (do_e) {
      e0t = 0;
      e0g = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec v0txx, v0tyx, v0tzx, v0tyy, v0tzy, v0tzz; // etors
   vbuf_prec v0gxx, v0gyx, v0gzx, v0gyy, v0gzy, v0gzz; // egeom
   if CONSTEXPR (do_v) {
      v0txx = 0, v0tyx = 0, v0tzx = 0, v0tyy = 0, v0tzy = 0, v0tzz = 0;
      v0gxx = 0, v0gyx = 0, v0gzx = 0, v0gyy = 0, v0gzy = 0, v0gzz = 0;
   }


   // etors
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
   if CONSTEXPR (do_e and rc_a) {
      if (ntors > 0)
         atomic_add(e0t, et, ithread);
   }
   if CONSTEXPR (do_v and rc_a) {
      if (ntors > 0)
         atomic_add(v0txx, v0tyx, v0tzx, v0tyy, v0tzy, v0tzz, vir_et, ithread);
   }


   // egeom
   for (int i = ithread; i < ngfix; i += stride) {
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_geom<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

                   degx, degy, degz,

                   i, igfix, gfix,

                   x, y, z, mass, molec, igrp, kgrp, grpmass,
                   TINKER_IMAGE_ARGS);
      if CONSTEXPR (do_e) {
         e0g += cvt_to<ebuf_prec>(e);
      }
      if CONSTEXPR (do_v) {
         v0gxx += cvt_to<vbuf_prec>(vxx);
         v0gyx += cvt_to<vbuf_prec>(vyx);
         v0gzx += cvt_to<vbuf_prec>(vzx);
         v0gyy += cvt_to<vbuf_prec>(vyy);
         v0gzy += cvt_to<vbuf_prec>(vzy);
         v0gzz += cvt_to<vbuf_prec>(vzz);
      }
   }
   if CONSTEXPR (do_e and rc_a) {
      if (ngfix > 0)
         atomic_add(e0g, eg, ithread);
   }
   if CONSTEXPR (do_v and rc_a) {
      if (ngfix > 0)
         atomic_add(v0gxx, v0gyx, v0gzx, v0gyy, v0gzy, v0gzz, vir_eg, ithread);
   }


   // total energy and virial
   if CONSTEXPR (do_e and not rc_a) {
      ebuf_prec etl = 0;
      etl += e0t; // etors
      etl += e0g; // egeom
      atomic_add(etl, ebuf, ithread);
   }
   if CONSTEXPR (do_v and not rc_a) {
      vbuf_prec vtlxx = 0, vtlyx = 0, vtlzx = 0;
      vbuf_prec vtlyy = 0, vtlzy = 0, vtlzz = 0;
      // etors
      vtlxx += v0txx, vtlyx += v0tyx, vtlzx += v0tzx;
      vtlyy += v0tyy, vtlzy += v0tzy, vtlzz += v0tzz;
      // egeom
      vtlxx += v0gxx, vtlyx += v0gyx, vtlzx += v0gzx;
      vtlyy += v0gyy, vtlzy += v0gzy, vtlzz += v0gzz;
      atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vbuf, ithread);
   }
}


void evalence_cu2(int vers, bool flag_tors, bool flag_geom)
{
#define EVALENCE_ARGS                                                          \
   /* etors */ et, vir_et, detx, dety, detz, torsunit, flag_tors ? ntors : 0,  \
      itors, tors1, tors2, tors3, tors4, tors5, tors6, /* egeom */ eg, vir_eg, \
      degx, degy, degz, flag_geom ? ngfix : 0, igfix, gfix,                    \
      /* total */ eng_buf, vir_buf, /* other */ x, y, z, mass,                 \
      molecule.molecule, grp.igrp, grp.kgrp, grp.grpmass, TINKER_IMAGE_ARGS


   int ngrid = get_grid_size(BLOCK_DIM);
   if (rc_flag & calc::analyz) {
      if (vers == calc::v0)
         evalence_cu1<calc::V0, true>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
      else if (vers == calc::v1)
         evalence_cu1<calc::V1, true>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
      else if (vers == calc::v3)
         evalence_cu1<calc::V3, true>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
      else if (vers == calc::v4)
         evalence_cu1<calc::V4, true>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
      else if (vers == calc::v5)
         evalence_cu1<calc::V5, true>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
      else if (vers == calc::v6)
         evalence_cu1<calc::V6, true>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
   } else {
      if (vers == calc::v0)
         evalence_cu1<calc::V0, false>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
      else if (vers == calc::v1)
         evalence_cu1<calc::V1, false>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
      else if (vers == calc::v3)
         assert(false);
      else if (vers == calc::v4)
         evalence_cu1<calc::V4, false>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
      else if (vers == calc::v5)
         evalence_cu1<calc::V5, false>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
      else if (vers == calc::v6)
         evalence_cu1<calc::V6, false>
            <<<ngrid, BLOCK_DIM, 0, nonblk>>>(EVALENCE_ARGS);
   }
}


void evalence_cu(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   bool flag_tors = use_potent(torsion_term);
   bool flag_geom = use_potent(geom_term);


   size_t bsize = buffer_size();
   if (rc_a and flag_tors) {
      host_zero(energy_et, virial_et);
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, et);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_et);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, detx, dety, detz);
   }
   if (rc_a and flag_geom) {
      host_zero(energy_eg, virial_eg);
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, eg);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_eg);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, degx, degy, degz);
   }


   if (flag_tors or flag_geom) {
      evalence_cu2(vers, flag_tors, flag_geom);
   }


   if (rc_a and flag_tors) {
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
   if (rc_a and flag_geom) {
      if (do_e) {
         energy_eg = energy_reduce(eg);
         energy_valence += energy_eg;
      }
      if (do_v) {
         virial_reduce(virial_eg, vir_eg);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eg[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, degx, degy, degz);
   }
}
}
