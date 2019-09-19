#include "acc_add.h"
#include "e_pitors.h"
#include "md.h"

TINKER_NAMESPACE_BEGIN
template <int USE>
void epitors_tmpl() {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  sanity_check<USE>();

  auto* ept = ept_handle.e()->buffer();
  auto* vir_ept = ept_handle.vir()->buffer();
  auto bufsize = ept_handle.buffer_size();

  #pragma acc parallel num_gangs(bufsize)\
              deviceptr(x,y,z,gx,gy,gz,\
              ipit,kpit,\
              ept,vir_ept)
  #pragma acc loop gang independent
  for (int i = 0; i < npitors; ++i) {
    int offset = i & (bufsize - 1);
    const int ia = ipit[i][0];
    const int ib = ipit[i][1];
    const int ic = ipit[i][2];
    const int id = ipit[i][3];
    const int ie = ipit[i][4];
    const int ig = ipit[i][5];

    real xia = x[ia];
    real yia = y[ia];
    real zia = z[ia];
    real xib = x[ib];
    real yib = y[ib];
    real zib = z[ib];
    real xic = x[ic];
    real yic = y[ic];
    real zic = z[ic];
    real xid = x[id];
    real yid = y[id];
    real zid = z[id];
    real xie = x[ie];
    real yie = y[ie];
    real zie = z[ie];
    real xig = x[ig];
    real yig = y[ig];
    real zig = z[ig];
    real xad = xia - xid;
    real yad = yia - yid;
    real zad = zia - zid;
    real xbd = xib - xid;
    real ybd = yib - yid;
    real zbd = zib - zid;
    real xec = xie - xic;
    real yec = yie - yic;
    real zec = zie - zic;
    real xgc = xig - xic;
    real ygc = yig - yic;
    real zgc = zig - zic;

    real xip = yad * zbd - ybd * zad + xic;
    real yip = zad * xbd - zbd * xad + yic;
    real zip = xad * ybd - xbd * yad + zic;
    real xiq = yec * zgc - ygc * zec + xid;
    real yiq = zec * xgc - zgc * xec + yid;
    real ziq = xec * ygc - xgc * yec + zid;
    real xcp = xic - xip;
    real ycp = yic - yip;
    real zcp = zic - zip;
    real xdc = xid - xic;
    real ydc = yid - yic;
    real zdc = zid - zic;
    real xqd = xiq - xid;
    real yqd = yiq - yid;
    real zqd = ziq - zid;

    real xt = ycp * zdc - ydc * zcp;
    real yt = zcp * xdc - zdc * xcp;
    real zt = xcp * ydc - xdc * ycp;
    real xu = ydc * zqd - yqd * zdc;
    real yu = zdc * xqd - zqd * xdc;
    real zu = xdc * yqd - xqd * ydc;
    real xtu = yt * zu - yu * zt;
    real ytu = zt * xu - zu * xt;
    real ztu = xt * yu - xu * yt;
    real rt2 = xt * xt + yt * yt + zt * zt;
    real ru2 = xu * xu + yu * yu + zu * zu;
    real rtru = REAL_SQRT(rt2 * ru2);
    if (rtru != 0) {
      real rdc = REAL_SQRT(xdc * xdc + ydc * ydc + zdc * zdc);
      real cosine = (xt * xu + yt * yu + zt * zu) * REAL_RECIP(rtru);
      real sine = (xdc * xtu + ydc * ytu + zdc * ztu) * REAL_RECIP(rdc * rtru);

      // set the pi-system torsion parameters for this angle

      real v2 = kpit[i];
      real c2 = -1;
      real s2 = 0;

      // compute the multiple angle trigonometry and the phase terms

      real cosine2 = cosine * cosine - sine * sine;
      real sine2 = 2 * cosine * sine;
      real phi2 = 1 + (cosine2 * c2 + sine2 * s2);

      // calculate the pi-system torsion energy for this angle

      if_constexpr(do_e) {
        real e = ptorunit * v2 * phi2;
        atomic_add_value(e, ept, offset);
      }

      if_constexpr(do_g) {
        real dphi2 = 2 * (cosine2 * s2 - sine2 * c2);
        real dedphi = ptorunit * v2 * dphi2;

        // chain rule terms for first derivative components

        real xdp = xid - xip;
        real ydp = yid - yip;
        real zdp = zid - zip;
        real xqc = xiq - xic;
        real yqc = yiq - yic;
        real zqc = ziq - zic;
        real rt2rdc_inv = REAL_RECIP(rt2 * rdc);
        real ru2rdc_inv = REAL_RECIP(ru2 * rdc);
        real dedxt = dedphi * (yt * zdc - ydc * zt) * rt2rdc_inv;
        real dedyt = dedphi * (zt * xdc - zdc * xt) * rt2rdc_inv;
        real dedzt = dedphi * (xt * ydc - xdc * yt) * rt2rdc_inv;
        real dedxu = -dedphi * (yu * zdc - ydc * zu) * ru2rdc_inv;
        real dedyu = -dedphi * (zu * xdc - zdc * xu) * ru2rdc_inv;
        real dedzu = -dedphi * (xu * ydc - xdc * yu) * ru2rdc_inv;

        // compute first derivative components for pi-system angle

        real dedxip = zdc * dedyt - ydc * dedzt;
        real dedyip = xdc * dedzt - zdc * dedxt;
        real dedzip = ydc * dedxt - xdc * dedyt;
        real dedxic = ydp * dedzt - zdp * dedyt + zqd * dedyu - yqd * dedzu;
        real dedyic = zdp * dedxt - xdp * dedzt + xqd * dedzu - zqd * dedxu;
        real dedzic = xdp * dedyt - ydp * dedxt + yqd * dedxu - xqd * dedyu;
        real dedxid = zcp * dedyt - ycp * dedzt + yqc * dedzu - zqc * dedyu;
        real dedyid = xcp * dedzt - zcp * dedxt + zqc * dedxu - xqc * dedzu;
        real dedzid = ycp * dedxt - xcp * dedyt + xqc * dedyu - yqc * dedxu;
        real dedxiq = zdc * dedyu - ydc * dedzu;
        real dedyiq = xdc * dedzu - zdc * dedxu;
        real dedziq = ydc * dedxu - xdc * dedyu;

        // compute first derivative components for individual atoms

        real dedxia = ybd * dedzip - zbd * dedyip;
        real dedyia = zbd * dedxip - xbd * dedzip;
        real dedzia = xbd * dedyip - ybd * dedxip;
        real dedxib = zad * dedyip - yad * dedzip;
        real dedyib = xad * dedzip - zad * dedxip;
        real dedzib = yad * dedxip - xad * dedyip;
        real dedxie = ygc * dedziq - zgc * dedyiq;
        real dedyie = zgc * dedxiq - xgc * dedziq;
        real dedzie = xgc * dedyiq - ygc * dedxiq;
        real dedxig = zec * dedyiq - yec * dedziq;
        real dedyig = xec * dedziq - zec * dedxiq;
        real dedzig = yec * dedxiq - xec * dedyiq;

        dedxic += (dedxip - dedxie - dedxig);
        dedyic += (dedyip - dedyie - dedyig);
        dedzic += (dedzip - dedzie - dedzig);
        dedxid += (dedxiq - dedxia - dedxib);
        dedyid += (dedyiq - dedyia - dedyib);
        dedzid += (dedziq - dedzia - dedzib);

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
        #pragma acc atomic update
        gx[id] += dedxid;
        #pragma acc atomic update
        gy[id] += dedyid;
        #pragma acc atomic update
        gz[id] += dedzid;
        #pragma acc atomic update
        gx[ie] += dedxie;
        #pragma acc atomic update
        gy[ie] += dedyie;
        #pragma acc atomic update
        gz[ie] += dedzie;
        #pragma acc atomic update
        gx[ig] += dedxig;
        #pragma acc atomic update
        gy[ig] += dedyig;
        #pragma acc atomic update
        gz[ig] += dedzig;

        if_constexpr(do_v) {
          real vxterm = dedxid + dedxia + dedxib;
          real vyterm = dedyid + dedyia + dedyib;
          real vzterm = dedzid + dedzia + dedzib;
          real vxx = xdc * vxterm + xcp * dedxip - xqd * dedxiq;
          real vyx = ydc * vxterm + ycp * dedxip - yqd * dedxiq;
          real vzx = zdc * vxterm + zcp * dedxip - zqd * dedxiq;
          real vyy = ydc * vyterm + ycp * dedyip - yqd * dedyiq;
          real vzy = zdc * vyterm + zcp * dedyip - zqd * dedyiq;
          real vzz = zdc * vzterm + zcp * dedzip - zqd * dedziq;

          int offv = offset * 8;
          atomic_add_value(vxx, vir_ept, offv + 0);
          atomic_add_value(vyx, vir_ept, offv + 1);
          atomic_add_value(vzx, vir_ept, offv + 2);
          atomic_add_value(vyy, vir_ept, offv + 3);
          atomic_add_value(vzy, vir_ept, offv + 4);
          atomic_add_value(vzz, vir_ept, offv + 5);
        }
      }
    }
  } // end for (int i)
}

void epitors_acc_impl_(int vers) {
  if (vers == calc::v0 || vers == calc::v3)
    epitors_tmpl<calc::v0>();
  else if (vers == calc::v1)
    epitors_tmpl<calc::v1>();
  else if (vers == calc::v4)
    epitors_tmpl<calc::v4>();
  else if (vers == calc::v5)
    epitors_tmpl<calc::v5>();
  else if (vers == calc::v6)
    epitors_tmpl<calc::v6>();
}
TINKER_NAMESPACE_END
