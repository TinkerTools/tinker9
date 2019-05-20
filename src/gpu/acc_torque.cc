#include "acc_e.h"
#include "gpu/e_mpole.h"

#define SET_ARRAY3_(a, ia, ib)                                                 \
  a[0] = x[ia] - x[ib];                                                        \
  a[1] = y[ia] - y[ib];                                                        \
  a[2] = z[ia] - z[ib]
#define NORM_(a) REAL_SQRT(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
#define CROSS_(ans, u, v)                                                      \
  ans[0] = u[1] * v[2] - u[2] * v[1];                                          \
  ans[1] = u[2] * v[0] - u[0] * v[2];                                          \
  ans[2] = u[0] * v[1] - u[1] * v[0]
#define NORMALIZE_(a, _1_na)                                                   \
  a[0] *= _1_na;                                                               \
  a[1] *= _1_na;                                                               \
  a[2] *= _1_na
#define ADD_(ans, a, b)                                                        \
  ans[0] = a[0] + b[0];                                                        \
  ans[1] = a[1] + b[1];                                                        \
  ans[2] = a[2] + b[2]
#define DOT_(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int DO_V>
void torque_tmpl(real* gpu_vir) {
  #pragma acc data deviceptr(x,y,z,gx,gy,gz,\
                             zaxis,trqx,trqy,trqz,\
                             gpu_vir)
  #pragma acc parallel loop
  for (int i = 0; i < n; ++i) {
    real frcz[3] = {0, 0, 0};
    real frcx[3] = {0, 0, 0};
    real frcy[3] = {0, 0, 0};
    const int axetyp = zaxis[i].polaxe;
    if (axetyp == pole_none)
      continue;

    const real trq[3] = {trqx[i], trqy[i], trqz[i]};
    real* de[3] = {gx, gy, gz};

    // get the local frame type and the frame-defining atoms

    const int ia = zaxis[i].zaxis;
    const int ib = i;
    const int ic = zaxis[i].xaxis;
    const int id = zaxis[i].yaxis;

    // construct the three rotation axes for the local frame

    real u[3], v[3], w[3];
    real usiz, vsiz, wsiz;
    real usiz1, vsiz1, wsiz1;
    SET_ARRAY3_(u, ia, ib);
    usiz = NORM_(u);
    usiz1 = REAL_RECIP(usiz);
    if (axetyp != pole_z_only) {
      SET_ARRAY3_(v, ic, ib);
      vsiz = NORM_(v);
      vsiz1 = REAL_RECIP(vsiz);
    } else {
      if (REAL_ABS(u[0]) > ((real)0.866) * usiz) {
        v[0] = 0;
        v[1] = 1;
      } else {
        v[0] = 1;
        v[1] = 0;
      }
      v[2] = 0;
      vsiz = 1;
      vsiz1 = 1;
    }
    if (axetyp == pole_z_bisect || axetyp == pole_3_fold) {
      SET_ARRAY3_(w, id, ib);
    } else {
      CROSS_(w, u, v);
    }
    wsiz = NORM_(w);
    wsiz1 = REAL_RECIP(wsiz);
    NORMALIZE_(u, usiz1);
    NORMALIZE_(v, vsiz1);
    NORMALIZE_(w, wsiz1);

    // build some additional axes needed for the Z-Bisect method

    real r[3], s[3];
    real rsiz, ssiz;
    real rsiz1, ssiz1;
    if (axetyp == pole_z_bisect) {
      ADD_(r, v, w);
      rsiz = NORM_(r);
      rsiz1 = REAL_RECIP(rsiz);
      CROSS_(s, u, r);
      ssiz = NORM_(s);
      ssiz1 = REAL_RECIP(ssiz);
      NORMALIZE_(r, rsiz1);
      NORMALIZE_(s, ssiz1);
    }

    // find the perpendicular and angle for each pair of axes

    real uv[3], uw[3], vw[3];
    real uvsiz, uwsiz, vwsiz;
    real uvsiz1, uwsiz1, vwsiz1;

    CROSS_(uv, v, u);
    uvsiz = NORM_(uv);
    uvsiz1 = REAL_RECIP(uvsiz);
    NORMALIZE_(uv, uvsiz1);

    CROSS_(uw, w, u);
    uwsiz = NORM_(uw);
    uwsiz1 = REAL_RECIP(uwsiz);
    NORMALIZE_(uw, uwsiz1);

    CROSS_(vw, w, v);
    vwsiz = NORM_(vw);
    vwsiz1 = REAL_RECIP(vwsiz);
    NORMALIZE_(vw, vwsiz1);

    real ur[3], us[3], vs[3], ws[3];
    real ursiz, ussiz, vssiz, wssiz;
    real ursiz1, ussiz1, vssiz1, wssiz1;
    if (axetyp == pole_z_bisect) {
      CROSS_(ur, r, u);
      ursiz = NORM_(ur);
      ursiz1 = REAL_RECIP(ursiz);
      NORMALIZE_(ur, ursiz1);

      CROSS_(us, s, u);
      ussiz = NORM_(us);
      ussiz1 = REAL_RECIP(ussiz);
      NORMALIZE_(us, ussiz1);

      CROSS_(vs, s, v);
      vssiz = NORM_(vs);
      vssiz1 = REAL_RECIP(vssiz);
      NORMALIZE_(vs, vssiz1);

      CROSS_(ws, s, w);
      wssiz = NORM_(ws);
      wssiz1 = REAL_RECIP(wssiz);
      NORMALIZE_(ws, wssiz1);
    }

    // get sine and cosine of angles between the rotation axes

    real uvcos = DOT_(u, v);
    real uvsin = REAL_SQRT(1 - uvcos * uvcos);
    real uvsin1 = REAL_RECIP(uvsin);
    real uwcos = DOT_(u, w);
    real uwsin = REAL_SQRT(1 - uwcos * uwcos);
    real vwcos = DOT_(v, w);
    real vwsin = REAL_SQRT(1 - vwcos * vwcos);

    real urcos, ursin, ursin1;
    real vscos, vssin;
    real wscos, wssin;
    if (axetyp == pole_z_bisect) {
      urcos = DOT_(u, r);
      ursin = REAL_SQRT(1 - urcos * urcos);
      ursin1 = REAL_RECIP(ursin);
      vscos = DOT_(v, s);
      vssin = REAL_SQRT(1 - vscos * vscos);
      wscos = DOT_(w, s);
      wssin = REAL_SQRT(1 - wscos * wscos);
    }

    // compute the projection of v and w onto the ru-plane

    real t1[3], t2[3];
    real t1siz, t2siz;
    real t1siz1, t2siz1;
    real ut1cos, ut1sin, ut2cos, ut2sin;
    if (axetyp == pole_z_bisect) {
      t1[0] = v[0] - s[0] * vscos;
      t1[1] = v[1] - s[1] * vscos;
      t1[2] = v[2] - s[2] * vscos;
      t1siz = NORM_(t1);
      t1siz1 = REAL_RECIP(t1siz);
      NORMALIZE_(t1, t1siz1);
      t2[0] = w[0] - s[0] * wscos;
      t2[1] = w[1] - s[1] * wscos;
      t2[2] = w[2] - s[2] * wscos;
      t2siz = NORM_(t2);
      t2siz1 = REAL_RECIP(t2siz);
      NORMALIZE_(t2, t2siz1);

      ut1cos = DOT_(u, t1);
      ut1sin = REAL_SQRT(1 - ut1cos * ut1cos);
      ut2cos = DOT_(u, t2);
      ut2sin = REAL_SQRT(1 - ut2cos * ut2cos);
    }

    // negative of dot product of torque with unit vectors gives result of
    // infinitesimal rotation along these vectors

    real dphidu = -DOT_(trq, u);
    real dphidv = -DOT_(trq, v);
    real dphidw = -DOT_(trq, w);
    real dphidr, dphids;
    if (axetyp == pole_z_bisect) {
      dphidr = -DOT_(trq, r);
      dphids = -DOT_(trq, s);
    }

    // force distribution
    if (axetyp == pole_z_only) {
      for (int j = 0; j < 3; ++j) {
        real du = uv[j] * dphidv * usiz1 * uvsin1 + uw[j] * dphidw * usiz1;
        #pragma acc atomic update
        de[j][ia] += du;
        #pragma acc atomic update
        de[j][ib] -= du;
        frcz[j] += du;
      }
    } else if (axetyp == pole_z_then_x) {
      for (int j = 0; j < 3; ++j) {
        real du = uv[j] * dphidv * usiz1 * uvsin1 + uw[j] * dphidw * usiz1;
        real dv = -uv[j] * dphidu * vsiz1 * uvsin1;
        #pragma acc atomic update
        de[j][ia] += du;
        #pragma acc atomic update
        de[j][ic] += dv;
        #pragma acc atomic update
        de[j][ib] -= (du + dv);
        frcz[j] += du;
        frcx[j] += dv;
      }
    } else if (axetyp == pole_bisector) {
      for (int j = 0; j < 3; ++j) {
        real du =
            uv[j] * dphidv * usiz1 * uvsin1 + 0.5f * uw[j] * dphidw * usiz1;
        real dv =
            -uv[j] * dphidu * vsiz1 * uvsin1 + 0.5f * vw[j] * dphidw * vsiz1;
        #pragma acc atomic update
        de[j][ia] += du;
        #pragma acc atomic update
        de[j][ic] += dv;
        #pragma acc atomic update
        de[j][ib] -= (du + dv);
        frcz[j] += du;
        frcx[j] += dv;
      }
    } else if (axetyp == pole_z_bisect) {
      for (int j = 0; j < 3; ++j) {
        real du = ur[j] * dphidr * usiz1 * ursin1 + us[j] * dphids * usiz1;
        real dv = (vssin * s[j] - vscos * t1[j]) * dphidu * vsiz1 *
            REAL_RECIP(ut1sin + ut2sin);
        real dw = (wssin * s[j] - wscos * t2[j]) * dphidu * wsiz1 *
            REAL_RECIP(ut1sin + ut2sin);
        #pragma acc atomic update
        de[j][ia] += du;
        #pragma acc atomic update
        de[j][ic] += dv;
        #pragma acc atomic update
        de[j][id] += dw;
        #pragma acc atomic update
        de[j][ib] -= (du + dv + dw);
        frcz[j] += du;
        frcx[j] += dv;
        frcy[j] += dw;
      }
    } else if (axetyp == pole_3_fold) {
      real p[3], psiz1;
      p[0] = u[0] + v[0] + w[0];
      p[1] = u[1] + v[1] + w[1];
      p[2] = u[2] + v[2] + w[2];
      psiz1 = REAL_RSQRT(DOT_(p, p));
      NORMALIZE_(p, psiz1);

      real wpcos;
      real upcos;
      real vpcos;
      wpcos = DOT_(w, p);
      upcos = DOT_(u, p);
      vpcos = DOT_(v, p);

      real r[3], rsiz1;
      real dphidr;
      real del[3], delsiz1;
      real dphiddel;
      real eps[3];

      ADD_(r, u, v);
      rsiz1 = REAL_RSQRT(DOT_(r, r));
      NORMALIZE_(r, rsiz1);
      real rwcos = DOT_(r, w);
      real rwsin1 = REAL_RSQRT(1 - rwcos * rwcos);
      dphidr = -DOT_(trq, r);
      CROSS_(del, r, w);
      delsiz1 = REAL_RSQRT(DOT_(del, del));
      NORMALIZE_(del, delsiz1);
      dphiddel = -DOT_(trq, del);
      CROSS_(eps, del, w);
      for (int j = 0; j < 3; ++j) {
        real dw = del[j] * dphidr * wsiz1 * rwsin1 +
            eps[j] * dphiddel * wpcos * wsiz1 * psiz1;
        #pragma acc atomic update
        de[j][id] += dw;
        #pragma acc atomic update
        de[j][ib] -= dw;
        frcy[j] += dw;
      }

      ADD_(r, v, w);
      rsiz1 = REAL_RSQRT(DOT_(r, r));
      NORMALIZE_(r, rsiz1);
      real rucos = DOT_(r, u);
      real rusin1 = REAL_RSQRT(1 - rucos * rucos);
      dphidr = -DOT_(trq, r);
      CROSS_(del, r, u);
      delsiz1 = REAL_RSQRT(DOT_(del, del));
      NORMALIZE_(del, delsiz1);
      dphiddel = -DOT_(trq, del);
      CROSS_(eps, del, u);
      for (int j = 0; j < 3; ++j) {
        real du = del[j] * dphidr * usiz1 * rusin1 +
            eps[j] * dphiddel * upcos * usiz1 * psiz1;
        #pragma acc atomic update
        de[j][ia] += du;
        #pragma acc atomic update
        de[j][ib] -= du;
        frcz[j] += du;
      }

      ADD_(r, u, w);
      rsiz1 = REAL_RSQRT(DOT_(r, r));
      NORMALIZE_(r, rsiz1);
      real rvcos = DOT_(r, v);
      real rvsin1 = REAL_RSQRT(1 - rvcos * rvcos);
      dphidr = -DOT_(trq, r);
      CROSS_(del, r, v);
      delsiz1 = REAL_RSQRT(DOT_(del, del));
      NORMALIZE_(del, delsiz1);
      dphiddel = -DOT_(trq, del);
      CROSS_(eps, del, v);
      for (int j = 0; j < 3; ++j) {
        real dv = del[j] * dphidr * vsiz1 * rvsin1 +
            eps[j] * dphiddel * vpcos * vsiz1 * psiz1;
        #pragma acc atomic update
        de[j][ic] += dv;
        #pragma acc atomic update
        de[j][ib] -= dv;
        frcx[j] += dv;
      }
    }

    if_constexpr(DO_V) {
      int iaz = (ia == -1) ? i : ia;
      int iax = (ic == -1) ? i : ic;
      int iay = (id == -1) ? i : id;

      real xiz = x[iaz] - x[i];
      real yiz = y[iaz] - y[i];
      real ziz = z[iaz] - z[i];

      real xix = x[iax] - x[i];
      real yix = y[iax] - y[i];
      real zix = z[iax] - z[i];

      real xiy = x[iay] - x[i];
      real yiy = y[iay] - y[i];
      real ziy = z[iay] - z[i];

      real vxx = xix * frcx[0] + xiy * frcy[0] + xiz * frcz[0];
      real vxy = 0.5f *
          (yix * frcx[0] + yiy * frcy[0] + yiz * frcz[0] + xix * frcx[1] +
           xiy * frcy[1] + xiz * frcz[1]);
      real vxz = 0.5f *
          (zix * frcx[0] + ziy * frcy[0] + ziz * frcz[0] + xix * frcx[2] +
           xiy * frcy[2] + xiz * frcz[2]);
      real vyy = yix * frcx[1] + yiy * frcy[1] + yiz * frcz[1];
      real vyz = 0.5f *
          (zix * frcx[1] + ziy * frcy[1] + ziz * frcz[1] + yix * frcx[2] +
           yiy * frcy[2] + yiz * frcz[2]);
      real vzz = zix * frcx[2] + ziy * frcy[2] + ziz * frcz[2];

      #pragma acc atomic update
      gpu_vir[_xx] += vxx;
      #pragma acc atomic update
      gpu_vir[_yx] += vxy;
      #pragma acc atomic update
      gpu_vir[_zx] += vxz;
      #pragma acc atomic update
      gpu_vir[_xy] += vxy;
      #pragma acc atomic update
      gpu_vir[_yy] += vyy;
      #pragma acc atomic update
      gpu_vir[_zy] += vyz;
      #pragma acc atomic update
      gpu_vir[_xz] += vxz;
      #pragma acc atomic update
      gpu_vir[_yz] += vyz;
      #pragma acc atomic update
      gpu_vir[_zz] += vzz;
    } // end if_constexpr(DO_V)
  }   // end for (int i)
}

void torque0() { torque_tmpl<0>(nullptr); }

void torque1() {
  if (use_empole())
    torque_tmpl<1>(vir_em);
}
}
TINKER_NAMESPACE_END

#undef SET_ARRAY3_
#undef NORM_
#undef CROSS_
#undef NORMALIZE_
#undef ADD_
#undef DOT_
