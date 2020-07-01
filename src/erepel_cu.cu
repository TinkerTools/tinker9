#include "add.h"
#include "erepel.h"
#include "glob.spatial.h"
#include "image.h"
#include "launch.h"
#include "md.h"
#include "seq_damprep.h"
#include "seq_switch.h"
#include "switch.h"

namespace tinker {
template <class Ver>
__device__
void pair_repel(real r2, real rscale, real cut, real off, real3 dr, real sizi,
                real dmpi, real vali, real ci, real3 di, real qixx, real qixy,
                real qixz, real qiyy, real qiyz, real qizz, //
                real sizk, real dmpk, real valk, real ck, real3 dk, real qkxx,
                real qkxy, real qkxz, real qkyy, real qkyz, real qkzz,
                real3& restrict frci, real3& restrict frck,
                real3& restrict trqi, real3& restrict trqk, real& restrict etl,
                real& restrict vtlxx, real& restrict vtlxy,
                real& restrict vtlxz, real& restrict vtlyy,
                real& restrict vtlyz, real& restrict vtlzz)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;

   real cut2 = cut * cut;
   real r = REAL_SQRT(r2);
   real rInv = REAL_RECIP(r);
   real rr2 = rInv * rInv;
   real rr1 = rInv;
   real rr3 = rr1 * rr2;
   real rr5 = 3 * rr3 * rr2;
   real rr7 = 5 * rr5 * rr2;
   real rr9 = 7 * rr7 * rr2;
   real rr11;
   real e;

   if CONSTEXPR (do_g) {
      rr11 = 9 * rr9 * rr2;
   }

   // get damping coefficients for the Pauli repulsion energy

   real dmpik[6];

   if CONSTEXPR (!do_g) {
      damp_rep<9>(dmpik, r, rr1, r2, rr3, rr5, rr7, rr9, rr11, dmpi, dmpk);
   } else {
      damp_rep<11>(dmpik, r, rr1, r2, rr3, rr5, rr7, rr9, rr11, dmpi, dmpk);
   }
   // dmpik(1) == dmpik[0]
   // dmpik(3) == dmpik[1]
   // dmpik(5) == dmpik[2]
   // dmpik(7) == dmpik[3]
   // dmpik(9) == dmpik[4]
   // dmpik(11) == dmpik[5]

   // same as in emplar_cu
   real3 frc, trq1, trq2;
   if CONSTEXPR (do_g) {
      frc = make_real3(0, 0, 0);
      trq1 = make_real3(0, 0, 0);
      trq2 = make_real3(0, 0, 0);
   }

   real dir = dot3(di, dr);
   real3 qi_dr = matvec(qixx, qixy, qixz, qiyy, qiyz, qizz, dr);
   real qir = dot3(dr, qi_dr);
   real dkr = dot3(dk, dr);
   real3 qk_dr = matvec(qkxx, qkxy, qkxz, qkyy, qkyz, qkzz, dr);
   real qkr = dot3(dr, qk_dr);

   real dik = dot3(di, dk);
   real qik = dot3(qi_dr, qk_dr);
   real diqk = dot3(di, qk_dr);
   real dkqi = dot3(dk, qi_dr);
   real qiqk = dot3(qixx, qiyy, qizz, qkxx, qkyy, qkzz) +
      2 * dot3(qixy, qixz, qiyz, qkxy, qkxz, qkyz);
   //

   real term1 = vali * valk;
   real term2 = valk * dir - vali * dkr + dik;
   real term3 = vali * qkr + valk * qir - dir * dkr + 2 * (dkqi - diqk + qiqk);
   real term4 = dir * qkr - dkr * qir - 4 * qik;
   real term5 = qir * qkr;

   real sizik = sizi * sizk * rscale;
   real eterm = term1 * dmpik[0] + term2 * dmpik[1] + term3 * dmpik[2] +
      term4 * dmpik[3] + term5 * dmpik[4];

   // energy
   if CONSTEXPR (do_e) {
      e = sizik * eterm * rInv;
   }

   // gradient
   real de = term1 * dmpik[1] + term2 * dmpik[2] + term3 * dmpik[3] +
      term4 * dmpik[4] + term5 * dmpik[5];
   term1 = -valk * dmpik[1] + dkr * dmpik[2] - qkr * dmpik[3];
   term2 = vali * dmpik[1] + dir * dmpik[2] + qir * dmpik[3];
   term3 = 2 * dmpik[2];
   term4 = 2 * (-valk * dmpik[2] + dkr * dmpik[3] - qkr * dmpik[4]);
   term5 = 2 * (-vali * dmpik[2] - dir * dmpik[3] - qir * dmpik[4]);
   real term6 = 4 * dmpik[3];

   // few intermediates, same as in emplar
   real qix = qi_dr.x;
   real qiy = qi_dr.y;
   real qiz = qi_dr.z;
   real qkx = qk_dr.x;
   real qky = qk_dr.y;
   real qkz = qk_dr.z;

   real qixk = qixx * qkx + qixy * qky + qixz * qkz; // |Qi Qk r>
   real qiyk = qixy * qkx + qiyy * qky + qiyz * qkz;
   real qizk = qixz * qkx + qiyz * qky + qizz * qkz;
   real qkxi = qkxx * qix + qkxy * qiy + qkxz * qiz; // |Qk Qi r>
   real qkyi = qkxy * qix + qkyy * qiy + qkyz * qiz;
   real qkzi = qkxz * qix + qkyz * qiy + qkzz * qiz;

   real diqkx = di.x * qkxx + di.y * qkxy + di.z * qkxz; // |Qk Di>
   real diqky = di.x * qkxy + di.y * qkyy + di.z * qkyz;
   real diqkz = di.x * qkxz + di.y * qkyz + di.z * qkzz;
   real dkqix = dk.x * qixx + dk.y * qixy + dk.z * qixz; // |Qi Dk>
   real dkqiy = dk.x * qixy + dk.y * qiyy + dk.z * qiyz;
   real dkqiz = dk.x * qixz + dk.y * qiyz + dk.z * qizz;

   //
   real3 frc0;
   frc0.x = de * dr.x + term1 * di.x + term2 * dk.x + term3 * (diqkx - dkqix) +
      term4 * qix + term5 * qkx + term6 * (qixk + qkxi);
   frc0.y = de * dr.y + term1 * di.y + term2 * dk.y + term3 * (diqky - dkqiy) +
      term4 * qiy + term5 * qky + term6 * (qiyk + qkyi);
   frc0.z = de * dr.z + term1 * di.z + term2 * dk.z + term3 * (diqkz - dkqiz) +
      term4 * qiz + term5 * qkz + term6 * (qizk + qkzi);

   // frc += frc0 * rr1 + eterm * rr3 * dr;
   // frc *= sizik;

   frc.x = frc0.x * rr1 + eterm * rr3 * dr.x;
   frc.y = frc0.y * rr1 + eterm * rr3 * dr.y;
   frc.z = frc0.z * rr1 + eterm * rr3 * dr.z;

   frc.x = frc.x * sizik;
   frc.y = frc.y * sizik;
   frc.z = frc.z * sizik;

   // torque
   real dirx = di.y * dr.z - di.z * dr.y; // Di x r
   real diry = di.z * dr.x - di.x * dr.z;
   real dirz = di.x * dr.y - di.y * dr.x;
   real dkrx = dk.y * dr.z - dk.z * dr.y; // Dk x r
   real dkry = dk.z * dr.x - dk.x * dr.z;
   real dkrz = dk.x * dr.y - dk.y * dr.x;
   real dikx = di.y * dk.z - di.z * dk.y; // Di x Dk
   real diky = di.z * dk.x - di.x * dk.z;
   real dikz = di.x * dk.y - di.y * dk.x;

   real qirx = qiz * dr.y - qiy * dr.z; // r x (Qi r)
   real qiry = qix * dr.z - qiz * dr.x;
   real qirz = qiy * dr.x - qix * dr.y;
   real qkrx = qkz * dr.y - qky * dr.z; // r x (Qk r)
   real qkry = qkx * dr.z - qkz * dr.x;
   real qkrz = qky * dr.x - qkx * dr.y;
   real qikx = qky * qiz - qkz * qiy; // (Qk r) x (Qi r)
   real qiky = qkz * qix - qkx * qiz;
   real qikz = qkx * qiy - qky * qix;

   real qikrx = qizk * dr.y - qiyk * dr.z;
   real qikry = qixk * dr.z - qizk * dr.x;
   real qikrz = qiyk * dr.x - qixk * dr.y;
   real qkirx = qkzi * dr.y - qkyi * dr.z;
   real qkiry = qkxi * dr.z - qkzi * dr.x;
   real qkirz = qkyi * dr.x - qkxi * dr.y;

   real diqkrx = diqkz * dr.y - diqky * dr.z;
   real diqkry = diqkx * dr.z - diqkz * dr.x;
   real diqkrz = diqky * dr.x - diqkx * dr.y;
   real dkqirx = dkqiz * dr.y - dkqiy * dr.z;
   real dkqiry = dkqix * dr.z - dkqiz * dr.x;
   real dkqirz = dkqiy * dr.x - dkqix * dr.y;

   real dqikx = di.y * qkz - di.z * qky + dk.y * qiz - dk.z * qiy -
      2 *
         (qixy * qkxz + qiyy * qkyz + qiyz * qkzz - qixz * qkxy - qiyz * qkyy -
          qizz * qkyz);
   real dqiky = di.z * qkx - di.x * qkz + dk.z * qix - dk.x * qiz -
      2 *
         (qixz * qkxx + qiyz * qkxy + qizz * qkxz - qixx * qkxz - qixy * qkyz -
          qixz * qkzz);
   real dqikz = di.x * qky - di.y * qkx + dk.x * qiy - dk.y * qix -
      2 *
         (qixx * qkxy + qixy * qkyy + qixz * qkyz - qixy * qkxx - qiyy * qkxy -
          qiyz * qkxz);

   real3 trq0;
   trq0.x = -dmpik[1] * dikx + term1 * dirx + term3 * (dqikx + dkqirx) -
      term4 * qirx - term6 * (qikrx + qikx);
   trq0.y = -dmpik[1] * diky + term1 * diry + term3 * (dqiky + dkqiry) -
      term4 * qiry - term6 * (qikry + qiky);
   trq0.z = -dmpik[1] * dikz + term1 * dirz + term3 * (dqikz + dkqirz) -
      term4 * qirz - term6 * (qikrz + qikz);
   trq1 += sizik * rr1 * trq0;

   trq0.x = dmpik[1] * dikx + term2 * dkrx - term3 * (dqikx + diqkrx) -
      term5 * qkrx - term6 * (qkirx - qikx);
   trq0.y = dmpik[1] * diky + term2 * dkry - term3 * (dqiky + diqkry) -
      term5 * qkry - term6 * (qkiry - qiky);
   trq0.z = dmpik[1] * dikz + term2 * dkrz - term3 * (dqikz + diqkrz) -
      term5 * qkrz - term6 * (qkirz - qikz);
   trq2 += sizik * rr1 * trq0;

   // no scaling based on group membership ?
   if (r2 > cut2) {
      real taper, dtaper;
      real cut = REAL_SQRT(cut2);
      switch_taper5<do_g>(r, cut, off, taper, dtaper);
      if CONSTEXPR (do_g)
         dtaper = dtaper * e * rr1;
      frc = frc * taper - dtaper * dr;
      trq1 *= taper;
      trq2 *= taper;
      if CONSTEXPR (do_e)
         e *= taper;
   }

   // save the results


   if CONSTEXPR (do_e) {
      etl += e;
   }

   if CONSTEXPR (do_g) {
      frci += frc;
      frck -= frc;
      trqi += trq1;
      trqk += trq2;
   }

   if CONSTEXPR (do_v) {
      vtlxx -= dr.x * frc.x;
      vtlxy -= 0.5f * (dr.y * frc.x + dr.x * frc.y);
      vtlxz -= 0.5f * (dr.z * frc.x + dr.x * frc.z);
      vtlyy -= dr.y * frc.y;
      vtlyz -= 0.5f * (dr.z * frc.y + dr.y * frc.z);
      vtlzz -= dr.z * frc.z;
   }
}


//====================================================================//


#define HIPPO_REPEL_PARA                                                       \
   size_t bufsize, count_buffer restrict nrep, energy_buffer restrict er,      \
      virial_buffer restrict vir_er, grad_prec *restrict derx,                 \
      grad_prec *restrict dery, grad_prec *restrict derz, real *restrict trqx, \
      real *restrict trqy, real *restrict trqz, TINKER_IMAGE_PARAMS,           \
      real *restrict sizpr, real *restrict dmppr, real *restrict elepr,        \
      real cut, real off, const real(*restrict rpole)[10]


template <class Ver>
__global__
void erepel_cu1(HIPPO_REPEL_PARA, int n,
                const Spatial::SortedAtom* restrict sorted, int niak,
                const int* restrict iak, const int* restrict lst)
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
      real3 pos;
      real3 frc, trq; // force and torque
      real c;         // charge
      real3 d;        // dipole
      real qxx, qxy, qxz, qyy, qyz, qzz;
      real siz, dmp, val;
   };
   __shared__ Data data[BLOCK_DIM];

   const real off2 = off * off;
   for (int iw = iwarp; iw < niak; iw += nwarp) {
      real etl;
      int ctl;
      if CONSTEXPR (do_a)
         ctl = 0;
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
         idat.frc = make_real3(0, 0, 0);
         idat.trq = make_real3(0, 0, 0);
      }
      int atomi = min(iak[iw] * WARP_SIZE + ilane, n - 1);
      idat.pos = make_real3(sorted[atomi].x, sorted[atomi].y, sorted[atomi].z);
      int i = sorted[atomi].unsorted;
      idat.c = rpole[i][mpl_pme_0];
      idat.d = make_real3(rpole[i][mpl_pme_x], rpole[i][mpl_pme_y],
                          rpole[i][mpl_pme_z]);
      idat.qxx = rpole[i][mpl_pme_xx];
      idat.qxy = rpole[i][mpl_pme_xy];
      idat.qxz = rpole[i][mpl_pme_xz];
      idat.qyy = rpole[i][mpl_pme_yy];
      idat.qyz = rpole[i][mpl_pme_yz];
      idat.qzz = rpole[i][mpl_pme_zz];
      idat.siz = sizpr[i];
      idat.dmp = dmppr[i];
      idat.val = elepr[i];


      if CONSTEXPR (do_g) {
         data[threadIdx.x].frc = make_real3(0, 0, 0);
         data[threadIdx.x].trq = make_real3(0, 0, 0);
      }
      int shatomk = lst[iw * WARP_SIZE + ilane];
      data[threadIdx.x].pos =
         make_real3(sorted[shatomk].x, sorted[shatomk].y, sorted[shatomk].z);
      int shk = sorted[shatomk].unsorted;
      data[threadIdx.x].c = rpole[shk][mpl_pme_0];
      data[threadIdx.x].d = make_real3(
         rpole[shk][mpl_pme_x], rpole[shk][mpl_pme_y], rpole[shk][mpl_pme_z]);
      data[threadIdx.x].qxx = rpole[shk][mpl_pme_xx];
      data[threadIdx.x].qxy = rpole[shk][mpl_pme_xy];
      data[threadIdx.x].qxz = rpole[shk][mpl_pme_xz];
      data[threadIdx.x].qyy = rpole[shk][mpl_pme_yy];
      data[threadIdx.x].qyz = rpole[shk][mpl_pme_yz];
      data[threadIdx.x].qzz = rpole[shk][mpl_pme_zz];
      data[threadIdx.x].siz = sizpr[shk];
      data[threadIdx.x].dmp = dmppr[shk];
      data[threadIdx.x].val = elepr[shk];


      for (int j = 0; j < WARP_SIZE; ++j) {
         int srclane = (ilane + j) & (WARP_SIZE - 1);
         int klane = srclane + threadIdx.x - ilane;
         int atomk = __shfl_sync(ALL_LANES, shatomk, srclane);
         real3 dr = data[klane].pos - idat.pos;


         real r2 = image2(dr.x, dr.y, dr.z);
         if (atomi < atomk && r2 <= off2) {

            pair_repel<Ver>(
               r2, 1, cut, off, dr, idat.siz, idat.dmp, idat.val, idat.c,
               idat.d, idat.qxx, idat.qxy, idat.qxz, idat.qyy, idat.qyz,
               idat.qzz, data[klane].siz, data[klane].dmp, data[klane].val,
               data[klane].c, data[klane].d, data[klane].qxx, data[klane].qxy,
               data[klane].qxz, data[klane].qyy, data[klane].qyz,
               data[klane].qzz, idat.frc, data[klane].frc, idat.trq,
               data[klane].trq, etl, vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz);


            if CONSTEXPR (do_a) {
               ctl += 1;
            }

         } // end if include
      }    // end for j loop
      if CONSTEXPR (do_a)
         atomic_add(ctl, nrep, offset);
      if CONSTEXPR (do_e)
         atomic_add(etl, er, offset);
      if CONSTEXPR (do_g) {
         atomic_add(idat.frc.x, &derx[i]);
         atomic_add(idat.frc.y, &dery[i]);
         atomic_add(idat.frc.z, &derz[i]);
         atomic_add(data[threadIdx.x].frc.x, &derx[shk]);
         atomic_add(data[threadIdx.x].frc.y, &dery[shk]);
         atomic_add(data[threadIdx.x].frc.z, &derz[shk]);
         atomic_add(idat.trq.x, &trqx[i]);
         atomic_add(idat.trq.y, &trqy[i]);
         atomic_add(idat.trq.z, &trqz[i]);
         atomic_add(data[threadIdx.x].trq.x, &trqx[shk]);
         atomic_add(data[threadIdx.x].trq.y, &trqy[shk]);
         atomic_add(data[threadIdx.x].trq.z, &trqz[shk]);
      }
      if CONSTEXPR (do_v)
         atomic_add(vtlxx, vtlxy, vtlxz, vtlyy, vtlyz, vtlzz, vir_er, offset);
   } // end for (iw)
}

template <class Ver>
__global__
void erepel_cu2(HIPPO_REPEL_PARA, const real* x, const real* y, const real* z,
                int nrepexclude, int (*restrict repexclude)[2],
                real* restrict repexclude_scale)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;

   const real off2 = off * off;
   for (int ii = threadIdx.x + blockIdx.x * blockDim.x; ii < nrepexclude;
        ii += blockDim.x * gridDim.x) {
      int offset = ii & (bufsize - 1);

      int i = repexclude[ii][0];
      int k = repexclude[ii][1];
      real rscale = repexclude_scale[ii];

      real xi = x[i];
      real yi = y[i];
      real zi = z[i];
      real ci = rpole[i][mpl_pme_0];
      real dix = rpole[i][mpl_pme_x];
      real diy = rpole[i][mpl_pme_y];
      real diz = rpole[i][mpl_pme_z];
      real qixx = rpole[i][mpl_pme_xx];
      real qixy = rpole[i][mpl_pme_xy];
      real qixz = rpole[i][mpl_pme_xz];
      real qiyy = rpole[i][mpl_pme_yy];
      real qiyz = rpole[i][mpl_pme_yz];
      real qizz = rpole[i][mpl_pme_zz];
      real sizi = sizpr[i];
      real dmpi = dmppr[i];
      real vali = elepr[i];


      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;
      real r2 = image2(xr, yr, zr);
      if (r2 <= off2) {
         real e;

         if CONSTEXPR (do_e) {
            e = 0;
         }
         real vxx, vxy, vxz, vyy, vyz, vzz;
         if CONSTEXPR (do_v) {
            vxx = 0;
            vxy = 0;
            vxz = 0;
            vyy = 0;
            vyz = 0;
            vzz = 0;
         }
         real3 frci, frck, trqi, trqk;
         if CONSTEXPR (do_g) {
            frci = make_real3(0, 0, 0);
            frck = make_real3(0, 0, 0);
            trqi = make_real3(0, 0, 0);
            trqk = make_real3(0, 0, 0);
         }
         real ck = rpole[k][mpl_pme_0];
         real dkx = rpole[k][mpl_pme_x];
         real dky = rpole[k][mpl_pme_y];
         real dkz = rpole[k][mpl_pme_z];
         real qkxx = rpole[k][mpl_pme_xx];
         real qkxy = rpole[k][mpl_pme_xy];
         real qkxz = rpole[k][mpl_pme_xz];
         real qkyy = rpole[k][mpl_pme_yy];
         real qkyz = rpole[k][mpl_pme_yz];
         real qkzz = rpole[k][mpl_pme_zz];
         real sizk = sizpr[k];
         real dmpk = dmppr[k];
         real valk = elepr[k];


         pair_repel<Ver>(r2, rscale, cut, off, make_real3(xr, yr, zr), sizi,
                         dmpi, vali, ci, make_real3(dix, diy, diz), qixx, qixy,
                         qixz, qiyy, qiyz, qizz, sizk, dmpk, valk, ck,
                         make_real3(dkx, dky, dkz), qkxx, qkxy, qkxz, qkyy,
                         qkyz, qkzz, frci, frck, trqi, trqk, e, vxx, vxy, vxz,
                         vyy, vyz, vzz);

         if CONSTEXPR (do_a)
            if (rscale == -1 && e != 0)
               atomic_add(-1, nrep, offset);
         if CONSTEXPR (do_e)
            atomic_add(e, er, offset);
         if CONSTEXPR (do_g) {
            atomic_add(frci.x, &derx[i]);
            atomic_add(frci.y, &dery[i]);
            atomic_add(frci.z, &derz[i]);
            atomic_add(frck.x, &derx[k]);
            atomic_add(frck.y, &dery[k]);
            atomic_add(frck.z, &derz[k]);
            atomic_add(trqi.x, &trqx[i]);
            atomic_add(trqi.y, &trqy[i]);
            atomic_add(trqi.z, &trqz[i]);
            atomic_add(trqk.x, &trqx[k]);
            atomic_add(trqk.y, &trqy[k]);
            atomic_add(trqk.z, &trqz[k]);
         }
         if CONSTEXPR (do_v)
            atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, vir_er, offset);
      }
   }
}

template <class Ver>
void erepel_cu3()
{
   const auto& st = *mspatial_unit;
   real cut = switch_cut(switch_repuls);
   real off = switch_off(switch_repuls);

   auto bufsize = buffer_size();

   auto ker1 = erepel_cu1<Ver>;
   if (st.niak > 0)
      launch_k1s(nonblk, WARP_SIZE * st.niak, ker1, //
                 bufsize, nrep, er, vir_er, derx, dery, derz, trqx, trqy, trqz,
                 TINKER_IMAGE_ARGS, sizpr, dmppr, elepr, cut, off, rpole, n,
                 st.sorted, st.niak, st.iak, st.lst);


   auto ker2 = erepel_cu2<Ver>;
   if (nrepexclude > 0)
      launch_k1s(nonblk, nrepexclude, ker2, //
                 bufsize, nrep, er, vir_er, derx, dery, derz, trqx, trqy, trqz,
                 TINKER_IMAGE_ARGS, sizpr, dmppr, elepr, cut, off, rpole, x, y,
                 z, nrepexclude, repexclude, repexclude_scale);
}


void erepel_cu(int vers)
{
   if (vers == calc::v0)
      erepel_cu3<calc::V0>();
   else if (vers == calc::v1)
      erepel_cu3<calc::V1>();
   else if (vers == calc::v3)
      erepel_cu3<calc::V3>();
   else if (vers == calc::v4)
      erepel_cu3<calc::V4>();
   else if (vers == calc::v5)
      erepel_cu3<calc::V5>();
   else if (vers == calc::v6)
      erepel_cu3<calc::V6>();
}
}