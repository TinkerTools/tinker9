#include "acc_seq.h"
#include "e_angle.h"
#include "gpu/e_bond.h"
#include "gpu/e_strbnd.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
template <int USE>
void estrbnd_tmpl() {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  sanity_check<USE>();

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,gx,gy,gz,\
              isb,sbk,iang,anat,bl,\
              eba,vir_eba)
  for (int istrbnd = 0; istrbnd < nstrbnd; ++istrbnd) {
    int i = isb[istrbnd][0];
    int j = isb[istrbnd][1];
    int k = isb[istrbnd][2];
    int ia = iang[i][0];
    int ib = iang[i][1];
    int ic = iang[i][2];
    real force1 = sbk[istrbnd][0];
    real force2 = sbk[istrbnd][1];

    real xia = x[ia];
    real yia = y[ia];
    real zia = z[ia];
    real xib = x[ib];
    real yib = y[ib];
    real zib = z[ib];
    real xic = x[ic];
    real yic = y[ic];
    real zic = z[ic];

    real xab = xia - xib;
    real yab = yia - yib;
    real zab = zia - zib;
    real xcb = xic - xib;
    real ycb = yic - yib;
    real zcb = zic - zib;

    real rab2 = xab * xab + yab * yab + zab * zab;
    real rcb2 = xcb * xcb + ycb * ycb + zcb * zcb;

    if (REAL_MIN(rab2, rcb2) != 0) {
      real rab = REAL_SQRT(rab2);
      real rcb = REAL_SQRT(rcb2);
      real xp = ycb * zab - zcb * yab;
      real yp = zcb * xab - xcb * zab;
      real zp = xcb * yab - ycb * xab;
      real rp = REAL_SQRT(xp * xp + yp * yp + zp * zp);
      rp = REAL_MAX(rp, (real)(0.0001));
      real dot = xab * xcb + yab * ycb + zab * zcb;
      real cosine = dot * REAL_RECIP(rab * rcb);
      cosine = REAL_MIN(1, REAL_MAX(-1, cosine));
      real angle = radian * REAL_ACOS(cosine);

      real dt = angle - anat[i];
      real term1 = -radian * REAL_RECIP(rab2 * rp);
      real term2 = radian * REAL_RECIP(rcb2 * rp);
      real ddtdxia = term1 * (yab * zp - zab * yp);
      real ddtdyia = term1 * (zab * xp - xab * zp);
      real ddtdzia = term1 * (xab * yp - yab * xp);
      real ddtdxic = term2 * (ycb * zp - zcb * yp);
      real ddtdyic = term2 * (zcb * xp - xcb * zp);
      real ddtdzic = term2 * (xcb * yp - ycb * xp);

      real dr1 = rab - bl[j];
      term1 = REAL_RECIP(rab);
      real dr2 = rcb - bl[k];
      term2 = REAL_RECIP(rcb);
      real ddrdxia = term1 * xab;
      real ddrdyia = term1 * yab;
      real ddrdzia = term1 * zab;
      real ddrdxic = term2 * xcb;
      real ddrdyic = term2 * ycb;
      real ddrdzic = term2 * zcb;

      term1 = stbnunit * force1;
      term2 = stbnunit * force2;
      real termr = term1 * dr1 + term2 * dr2;

      if_constexpr(do_e) {
        real e = termr * dt;
        #pragma acc atomic update
        *eba += e;
      }

      if_constexpr(do_g) {
        real term1t = term1 * dt;
        real term2t = term2 * dt;
        real dedxia = term1t * ddrdxia + termr * ddtdxia;
        real dedyia = term1t * ddrdyia + termr * ddtdyia;
        real dedzia = term1t * ddrdzia + termr * ddtdzia;
        real dedxic = term2t * ddrdxic + termr * ddtdxic;
        real dedyic = term2t * ddrdyic + termr * ddtdyic;
        real dedzic = term2t * ddrdzic + termr * ddtdzic;
        real dedxib = -dedxia - dedxic;
        real dedyib = -dedyia - dedyic;
        real dedzib = -dedzia - dedzic;

        #pragma acc atomic update
        gx[ia] += dedxia;
        #pragma acc atomic update
        gy[ia] += dedyia;
        #pragma acc atomic update
        gz[ia] += dedzia;
        #pragma acc atomic update
        gx[ib] += dedxib;
        #pragma acc atomic update
        gy[ib] += dedyib;
        #pragma acc atomic update
        gz[ib] += dedzib;
        #pragma acc atomic update
        gx[ic] += dedxic;
        #pragma acc atomic update
        gy[ic] += dedyic;
        #pragma acc atomic update
        gz[ic] += dedzic;

        if_constexpr(do_v) {
          real vxx = xab * dedxia + xcb * dedxic;
          real vyx = yab * dedxia + ycb * dedxic;
          real vzx = zab * dedxia + zcb * dedxic;
          real vyy = yab * dedyia + ycb * dedyic;
          real vzy = zab * dedyia + zcb * dedyic;
          real vzz = zab * dedzia + zcb * dedzic;

          #pragma acc atomic update
          vir_eba[0] += vxx;
          #pragma acc atomic update
          vir_eba[1] += vyx;
          #pragma acc atomic update
          vir_eba[2] += vzx;
          #pragma acc atomic update
          vir_eba[3] += vyx;
          #pragma acc atomic update
          vir_eba[4] += vyy;
          #pragma acc atomic update
          vir_eba[5] += vzy;
          #pragma acc atomic update
          vir_eba[6] += vzx;
          #pragma acc atomic update
          vir_eba[7] += vzy;
          #pragma acc atomic update
          vir_eba[8] += vzz;
        }
      }
    }
  }
}

void estrbnd_acc_impl_(int vers) {
  if (vers == calc::v0 || vers == calc::v3)
    estrbnd_tmpl<calc::v0>();
  else if (vers == calc::v1)
    estrbnd_tmpl<calc::v1>();
  else if (vers == calc::v4)
    estrbnd_tmpl<calc::v4>();
  else if (vers == calc::v5)
    estrbnd_tmpl<calc::v5>();
  else if (vers == calc::v6)
    estrbnd_tmpl<calc::v6>();
}
TINKER_NAMESPACE_END
