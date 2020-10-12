#include "add.h"
#include "energy.h"
#include "glob.group.h"
#include "glob.molecule.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "potent.h"
#include "seq_angle.h"
#include "seq_bond.h"
#include "seq_geom.h"
#include "seq_imptor.h"
#include "seq_opbend.h"
#include "seq_pitors.h"
#include "seq_strbnd.h"
#include "seq_torsion.h"
#include "seq_tortor.h"
#include "seq_urey.h"
#include "tool/gpu_card.h"
#include "tool/host_zero.h"


namespace tinker {
template <class Ver, bool rc_a>
__global__
void evalence_cu1(
   // ebond
   energy_buffer restrict eb, virial_buffer restrict vir_eb,
   grad_prec* restrict debx, grad_prec* restrict deby, grad_prec* restrict debz,

   ebond_t bndtyp, real bndunit, int nbond, const int (*restrict ibnd)[2],
   const real* restrict bl, const real* restrict bk, real cbnd, real qbnd,

   // eimptor
   energy_buffer restrict eit, virial_buffer restrict vir_eit,
   grad_prec* restrict deitx, grad_prec* restrict deity,
   grad_prec* restrict deitz,

   real itorunit, int nitors, const int (*restrict iitors)[4],
   const real (*restrict itors1)[4], const real (*restrict itors2)[4],
   const real (*restrict itors3)[4],

   // etors
   energy_buffer restrict et, virial_buffer restrict vir_et,
   grad_prec* restrict detx, grad_prec* restrict dety, grad_prec* restrict detz,

   real torsunit, int ntors, const int (*restrict itors)[4],
   const real (*restrict tors1)[4], const real (*restrict tors2)[4],
   const real (*restrict tors3)[4], const real (*restrict tors4)[4],
   const real (*restrict tors5)[4], const real (*restrict tors6)[4],

   // epitors
   energy_buffer restrict ept, virial_buffer restrict vir_ept,
   grad_prec* restrict deptx, grad_prec* restrict depty,
   grad_prec* restrict deptz,

   real ptorunit, int npitors, const int (*restrict ipit)[6],
   const real* restrict kpit,

   // etortor
   energy_buffer restrict ett, virial_buffer restrict vir_ett,
   grad_prec* restrict dettx, grad_prec* restrict detty,
   grad_prec* restrict dettz,

   real ttorunit, int ntortor, const int (*restrict itt)[3],
   const int (*restrict ibitor)[5], const int* restrict chkttor_ia_,

   const int* restrict tnx, const int* restrict tny,
   const real (*restrict ttx)[ktrtor::maxtgrd],
   const real (*restrict tty)[ktrtor::maxtgrd],
   const real (*restrict tbf)[ktrtor::maxtgrd2],
   const real (*restrict tbx)[ktrtor::maxtgrd2],
   const real (*restrict tby)[ktrtor::maxtgrd2],
   const real (*restrict tbxy)[ktrtor::maxtgrd2],

   // egeom
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
   constexpr bool do_v = Ver::v;
   const int ithread = threadIdx.x + blockIdx.x * blockDim.x;
   const int stride = blockDim.x * gridDim.x;


   using ebuf_prec = energy_buffer_traits::type;
   ebuf_prec e0b;  // ebond
   ebuf_prec e0it; // eimptor
   ebuf_prec e0t;  // etors
   ebuf_prec e0pt; // epitors
   ebuf_prec e0tt; // etortor
   ebuf_prec e0g;  // egeom
   if CONSTEXPR (do_e) {
      e0b = 0;
      e0it = 0;
      e0t = 0;
      e0pt = 0;
      e0tt = 0;
      e0g = 0;
   }
   using vbuf_prec = virial_buffer_traits::type;
   vbuf_prec v0bxx, v0byx, v0bzx, v0byy, v0bzy, v0bzz;       // ebond
   vbuf_prec v0itxx, v0ityx, v0itzx, v0ityy, v0itzy, v0itzz; // eimptor
   vbuf_prec v0txx, v0tyx, v0tzx, v0tyy, v0tzy, v0tzz;       // etors
   vbuf_prec v0ptxx, v0ptyx, v0ptzx, v0ptyy, v0ptzy, v0ptzz; // epitors
   vbuf_prec v0ttxx, v0ttyx, v0ttzx, v0ttyy, v0ttzy, v0ttzz; // etors
   vbuf_prec v0gxx, v0gyx, v0gzx, v0gyy, v0gzy, v0gzz;       // egeom
   if CONSTEXPR (do_v) {
      v0bxx = 0, v0byx = 0, v0bzx = 0, v0byy = 0, v0bzy = 0, v0bzz = 0;
      v0itxx = 0, v0ityx = 0, v0itzx = 0, v0ityy = 0, v0itzy = 0, v0itzz = 0;
      v0txx = 0, v0tyx = 0, v0tzx = 0, v0tyy = 0, v0tzy = 0, v0tzz = 0;
      v0ptxx = 0, v0ptyx = 0, v0ptzx = 0, v0ptyy = 0, v0ptzy = 0, v0ptzz = 0;
      v0ttxx = 0, v0ttyx = 0, v0ttzx = 0, v0ttyy = 0, v0ttzy = 0, v0ttzz = 0;
      v0gxx = 0, v0gyx = 0, v0gzx = 0, v0gyy = 0, v0gzy = 0, v0gzz = 0;
   }


   // ebond
   for (int i = ithread; i < nbond; i += stride) {
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_bond<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

                   debx, deby, debz,

                   bndtyp, bndunit, i, ibnd, bl, bk, cbnd, qbnd,

                   x, y, z);
      if CONSTEXPR (do_e) {
         e0b += cvt_to<ebuf_prec>(e);
      }
      if CONSTEXPR (do_v) {
         v0bxx += cvt_to<vbuf_prec>(vxx);
         v0byx += cvt_to<vbuf_prec>(vyx);
         v0bzx += cvt_to<vbuf_prec>(vzx);
         v0byy += cvt_to<vbuf_prec>(vyy);
         v0bzy += cvt_to<vbuf_prec>(vzy);
         v0bzz += cvt_to<vbuf_prec>(vzz);
      }
   }
   if CONSTEXPR (do_e and rc_a) {
      if (nbond > 0)
         atomic_add(e0b, eb, ithread);
   }
   if CONSTEXPR (do_v and rc_a) {
      if (nbond > 0)
         atomic_add(v0bxx, v0byx, v0bzx, v0byy, v0bzy, v0bzz, vir_eb, ithread);
   }


   // eimptor
   for (int i = ithread; i < nitors; i += stride) {
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_imptor<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

                     deitx, deity, deitz,

                     itorunit, i, iitors, itors1, itors2, itors3,

                     x, y, z);
      if CONSTEXPR (do_e) {
         e0it += cvt_to<ebuf_prec>(e);
      }
      if CONSTEXPR (do_v) {
         v0itxx += cvt_to<vbuf_prec>(vxx);
         v0ityx += cvt_to<vbuf_prec>(vyx);
         v0itzx += cvt_to<vbuf_prec>(vzx);
         v0ityy += cvt_to<vbuf_prec>(vyy);
         v0itzy += cvt_to<vbuf_prec>(vzy);
         v0itzz += cvt_to<vbuf_prec>(vzz);
      }
   }
   if CONSTEXPR (do_e and rc_a) {
      if (nitors > 0)
         atomic_add(e0it, eit, ithread);
   }
   if CONSTEXPR (do_v and rc_a) {
      if (nitors > 0)
         atomic_add(v0itxx, v0ityx, v0itzx, v0ityy, v0itzy, v0itzz, vir_eit,
                    ithread);
   }


   // etors
   for (int i = ithread; i < ntors; i += stride) {
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_tors<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

                   detx, dety, detz,

                   torsunit, i, itors,

                   tors1, tors2, tors3, tors4, tors5, tors6, x, y, z);
      if CONSTEXPR (do_e) {
         e0t += cvt_to<ebuf_prec>(e);
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


   // epitors
   for (int i = ithread; i < npitors; i += stride) {
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_pitors<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

                     deptx, depty, deptz,

                     ptorunit, i, ipit, kpit, x, y, z);
      if CONSTEXPR (do_e) {
         e0pt += cvt_to<ebuf_prec>(e);
      }
      if CONSTEXPR (do_v) {
         v0ptxx += cvt_to<vbuf_prec>(vxx);
         v0ptyx += cvt_to<vbuf_prec>(vyx);
         v0ptzx += cvt_to<vbuf_prec>(vzx);
         v0ptyy += cvt_to<vbuf_prec>(vyy);
         v0ptzy += cvt_to<vbuf_prec>(vzy);
         v0ptzz += cvt_to<vbuf_prec>(vzz);
      }
   }
   if CONSTEXPR (do_e and rc_a) {
      if (npitors > 0)
         atomic_add(e0pt, ept, ithread);
   }
   if CONSTEXPR (do_v and rc_a) {
      if (npitors > 0)
         atomic_add(v0ptxx, v0ptyx, v0ptzx, v0ptyy, v0ptzy, v0ptzz, vir_ept,
                    ithread);
   }


   // etortor
   for (int i = ithread; i < ntortor; i += stride) {
      real e, vxx, vyx, vzx, vyy, vzy, vzz;
      dk_tortor<Ver>(e, vxx, vyx, vzx, vyy, vzy, vzz,

                     dettx, detty, dettz,

                     ttorunit, i, itt, ibitor, chkttor_ia_,

                     tnx, tny, ttx, tty, tbf, tbx, tby, tbxy,

                     x, y, z);
      if CONSTEXPR (do_e) {
         e0tt += cvt_to<ebuf_prec>(e);
      }
      if CONSTEXPR (do_v) {
         v0ttxx += cvt_to<vbuf_prec>(vxx);
         v0ttyx += cvt_to<vbuf_prec>(vyx);
         v0ttzx += cvt_to<vbuf_prec>(vzx);
         v0ttyy += cvt_to<vbuf_prec>(vyy);
         v0ttzy += cvt_to<vbuf_prec>(vzy);
         v0ttzz += cvt_to<vbuf_prec>(vzz);
      }
   }
   if CONSTEXPR (do_e and rc_a) {
      if (ntortor > 0)
         atomic_add(e0tt, ett, ithread);
   }
   if CONSTEXPR (do_v and rc_a) {
      if (ntortor > 0)
         atomic_add(v0ttxx, v0ttyx, v0ttzx, v0ttyy, v0ttzy, v0ttzz, vir_ett,
                    ithread);
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
      etl += e0b;  // ebond
      etl += e0it; // eimptor
      etl += e0t;  // etors
      etl += e0pt; // epitors
      etl += e0tt; // etortor
      etl += e0g;  // egeom
      atomic_add(etl, ebuf, ithread);
   }
   if CONSTEXPR (do_v and not rc_a) {
      vbuf_prec vtlxx = 0, vtlyx = 0, vtlzx = 0;
      vbuf_prec vtlyy = 0, vtlzy = 0, vtlzz = 0;
      // ebond
      vtlxx += v0bxx, vtlyx += v0byx, vtlzx += v0bzx;
      vtlyy += v0byy, vtlzy += v0bzy, vtlzz += v0bzz;
      // eimptor
      vtlxx += v0itxx, vtlyx += v0ityx, vtlzx += v0itzx;
      vtlyy += v0ityy, vtlzy += v0itzy, vtlzz += v0itzz;
      // etors
      vtlxx += v0txx, vtlyx += v0tyx, vtlzx += v0tzx;
      vtlyy += v0tyy, vtlzy += v0tzy, vtlzz += v0tzz;
      // epitors
      vtlxx += v0ptxx, vtlyx += v0ptyx, vtlzx += v0ptzx;
      vtlyy += v0ptyy, vtlzy += v0ptzy, vtlzz += v0ptzz;
      // etortor
      vtlxx += v0ttxx, vtlyx += v0ttyx, vtlzx += v0ttzx;
      vtlyy += v0ttyy, vtlzy += v0ttzy, vtlzz += v0ttzz;
      // egeom
      vtlxx += v0gxx, vtlyx += v0gyx, vtlzx += v0gzx;
      vtlyy += v0gyy, vtlzy += v0gzy, vtlzz += v0gzz;
      atomic_add(vtlxx, vtlyx, vtlzx, vtlyy, vtlzy, vtlzz, vbuf, ithread);
   }
}


void evalence_cu2(int vers, bool flag_bond, bool flag_angle, bool flag_strbnd,
                  bool flag_urey, bool flag_opb, bool flag_imptor,
                  bool flag_tors, bool flag_pitors, bool flag_tortor,
                  bool flag_geom)
{
#define EVALENCE_ARGS                                                          \
   /* ebond */ eb, vir_eb, debx, deby, debz, bndtyp, bndunit,                  \
      flag_bond ? nbond : 0, ibnd, bl, bk, cbnd, qbnd, /* eimptor */ eit,      \
      vir_eit, deitx, deity, deitz, itorunit, flag_imptor ? nitors : 0,        \
      iitors, itors1, itors2, itors3, /* etors */ et, vir_et, detx, dety,      \
      detz, torsunit, flag_tors ? ntors : 0, itors, tors1, tors2, tors3,       \
      tors4, tors5, tors6, /* epitors */ ept, vir_ept, deptx, depty, deptz,    \
      ptorunit, flag_pitors ? npitors : 0, ipit, kpit, /* etortor */ ett,      \
      vir_ett, dettx, detty, dettz, ttorunit, flag_tortor ? ntortor : 0, itt,  \
      ibitor, chkttor_ia_, tnx, tny, ttx, tty, tbf, tbx, tby, tbxy,            \
      /* egeom */ eg, vir_eg, degx, degy, degz, flag_geom ? ngfix : 0, igfix,  \
      gfix, /* total */ eng_buf, vir_buf, /* other */ x, y, z, mass,           \
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
#undef EVALENCE_ARGS
}


void evalence_cu(int vers)
{
   bool rc_a = rc_flag & calc::analyz;
   bool do_e = vers & calc::energy;
   bool do_v = vers & calc::virial;
   bool do_g = vers & calc::grad;


   bool flag_bond = use_potent(bond_term);
   bool flag_angle = use_potent(angle_term);
   bool flag_strbnd = use_potent(strbnd_term);
   bool flag_urey = use_potent(urey_term);
   bool flag_opb = use_potent(opbend_term);
   bool flag_imptor = use_potent(imptors_term);
   bool flag_tors = use_potent(torsion_term);
   bool flag_pitors = use_potent(pitors_term);
   bool flag_tortor = use_potent(tortor_term);
   bool flag_geom = use_potent(geom_term);


   size_t bsize = buffer_size();

   if (rc_a and flag_bond) {
      host_zero(energy_eb, virial_eb);
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, eb);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_eb);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, debx, deby, debz);
   }
   if (rc_a and flag_imptor) {
      size_t bsize = buffer_size();
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, eit);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_eit);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, deitx, deity, deitz);
   }
   if (rc_a and flag_tors) {
      host_zero(energy_et, virial_et);
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, et);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_et);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, detx, dety, detz);
   }
   if (rc_a and flag_pitors) {
      host_zero(energy_ept, virial_ept);
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, ept);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_ept);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, deptx, depty, deptz);
   }
   if (rc_a and flag_tortor) {
      host_zero(energy_ett, virial_ett);
      if (do_e)
         darray::zero(PROCEED_NEW_Q, bsize, ett);
      if (do_v)
         darray::zero(PROCEED_NEW_Q, bsize, vir_ett);
      if (do_g)
         darray::zero(PROCEED_NEW_Q, n, dettx, detty, dettz);
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


   if (flag_bond or flag_angle or flag_strbnd or flag_urey or flag_opb or
       flag_imptor or flag_tors or flag_pitors or flag_tortor or flag_geom) {
      evalence_cu2(vers, flag_bond, flag_angle, flag_strbnd, flag_urey,
                   flag_opb, flag_imptor, flag_tors, flag_pitors, flag_tortor,
                   flag_geom);
   }


   if (rc_a and flag_bond) {
      if (do_e) {
         energy_eb = energy_reduce(eb);
         energy_valence += energy_eb;
      }
      if (do_v) {
         virial_reduce(virial_eb, vir_eb);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eb[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, debx, deby, debz);
   }
   if (rc_a and flag_imptor) {
      if (do_e) {
         energy_eit = energy_reduce(eit);
         energy_valence += energy_eit;
      }
      if (do_v) {
         virial_reduce(virial_eit, vir_eit);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_eit[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, deitx, deity, deitz);
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
   if (rc_a and flag_pitors) {
      if (do_e) {
         energy_ept = energy_reduce(ept);
         energy_valence += energy_ept;
      }
      if (do_v) {
         virial_reduce(virial_ept, vir_ept);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ept[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, deptx, depty, deptz);
   }
   if (rc_a and flag_tortor) {
      if (do_e) {
         energy_ett = energy_reduce(ett);
         energy_valence += energy_ett;
      }
      if (do_v) {
         virial_reduce(virial_ett, vir_ett);
         for (int iv = 0; iv < 9; ++iv)
            virial_valence[iv] += virial_ett[iv];
      }
      if (do_g)
         sum_gradient(gx, gy, gz, dettx, detty, dettz);
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
