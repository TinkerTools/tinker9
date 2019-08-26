#include "acc_add.h"
#include "elec.h"
#include "mathfunc.h"
#include "md.h"

#define ADD_(ans, a, b)                                                        \
  ans[0] = a[0] + b[0];                                                        \
  ans[1] = a[1] + b[1];                                                        \
  ans[2] = a[2] + b[2]
#define SUB_(a, ia, ib)                                                        \
  a[0] = x[ia] - x[ib];                                                        \
  a[1] = y[ia] - y[ib];                                                        \
  a[2] = z[ia] - z[ib]
#define DOT_(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
#define CROSS_(ans, u, v)                                                      \
  ans[0] = u[1] * v[2] - u[2] * v[1];                                          \
  ans[1] = u[2] * v[0] - u[0] * v[2];                                          \
  ans[2] = u[0] * v[1] - u[1] * v[0]
#define NORMAL_(a, _1_na)                                                      \
  a[0] *= _1_na;                                                               \
  a[1] *= _1_na;                                                               \
  a[2] *= _1_na

TINKER_NAMESPACE_BEGIN
template <int DO_V>
void torque_tmpl(Virial v_handle) {
  VirialBuffer::PointerType gpu_vir = nullptr;
  auto bufsize = VirialBuffer::estimate_size(n);
  if (v_handle >= 0) {
    gpu_vir = v_handle->buffer();
    bufsize = v_handle->size();
  }

  #pragma acc parallel loop gang num_gangs(bufsize) independent\
              deviceptr(x,y,z,gx,gy,gz,zaxis,trqx,trqy,trqz,gpu_vir)
  for (int i = 0; i < n; ++i) {
    const int axetyp = zaxis[i].polaxe;
    if (axetyp == pole_none)
      continue;

    const real trq[3] = {trqx[i], trqy[i], trqz[i]};
    real* de[3] = {gx, gy, gz};

    // zero out force components on local frame-defining atoms

    real frcz[3] = {0, 0, 0};
    real frcx[3] = {0, 0, 0};
    real frcy[3] = {0, 0, 0};

    // get the local frame type and the frame-defining atoms

    const int ia = zaxis[i].zaxis;
    const int ib = i;
    const int ic = zaxis[i].xaxis;
    const int id = INT_ABS(zaxis[i].yaxis) - 1;

    // construct the three rotation axes for the local frame

    real u[3], usiz1;
    SUB_(u, ia, ib);
    usiz1 = REAL_RSQRT(DOT_(u, u));

    real v[3], vsiz1;
    if (axetyp != pole_z_only) {
      SUB_(v, ic, ib);
      vsiz1 = REAL_RSQRT(DOT_(v, v));
    } else {
      if (REAL_ABS(u[0]) * usiz1 > ((real)0.866)) {
        v[0] = 0;
        v[1] = 1;
      } else {
        v[0] = 1;
        v[1] = 0;
      }
      v[2] = 0;
      vsiz1 = 1;
    }

    real w[3], wsiz1;
    if (axetyp == pole_z_bisect || axetyp == pole_3_fold) {
      SUB_(w, id, ib);
    } else {
      CROSS_(w, u, v);
    }
    wsiz1 = REAL_RSQRT(DOT_(w, w));

    NORMAL_(u, usiz1);
    NORMAL_(v, vsiz1);
    NORMAL_(w, wsiz1);

    // find the perpendicular and angle for each pair of axes

    real uv[3], uvsiz1;
    CROSS_(uv, v, u);
    uvsiz1 = REAL_RSQRT(DOT_(uv, uv));
    NORMAL_(uv, uvsiz1);

    real uw[3], uwsiz1;
    CROSS_(uw, w, u);
    uwsiz1 = REAL_RSQRT(DOT_(uw, uw));
    NORMAL_(uw, uwsiz1);

    real vw[3], vwsiz1;
    CROSS_(vw, w, v);
    vwsiz1 = REAL_RSQRT(DOT_(vw, vw));
    NORMAL_(vw, vwsiz1);

    // get sine and cosine of angles between the rotation axes

    real uvcos = DOT_(u, v);
    real uvsin1 = REAL_RSQRT(1 - uvcos * uvcos);

    // negative of dot product of torque with unit vectors gives result of
    // infinitesimal rotation along these vectors

    real dphidu = -DOT_(trq, u);
    real dphidv = -DOT_(trq, v);
    real dphidw = -DOT_(trq, w);

    // force distribution

    if (axetyp == pole_z_only) {
      for (int j = 0; j < 3; ++j) {
        real du = uv[j] * dphidv * usiz1 * uvsin1 + uw[j] * dphidw * usiz1;
        #pragma acc atomic update
        de[j][ia] += du;
        #pragma acc atomic update
        de[j][ib] -= du;
        if_constexpr(DO_V) frcz[j] += du;
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
        if_constexpr(DO_V) {
          frcz[j] += du;
          frcx[j] += dv;
        }
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
        if_constexpr(DO_V) {
          frcz[j] += du;
          frcx[j] += dv;
        }
      }
    } else if (axetyp == pole_z_bisect) {

      // build some additional axes needed for the Z-Bisect method

      real r[3], rsiz1;
      ADD_(r, v, w);
      rsiz1 = REAL_RSQRT(DOT_(r, r));

      real s[3], ssiz1;
      CROSS_(s, u, r);
      ssiz1 = REAL_RSQRT(DOT_(s, s));

      NORMAL_(r, rsiz1);
      NORMAL_(s, ssiz1);

      // find the perpendicular and angle for each pair of axes

      real ur[3], ursiz1;
      CROSS_(ur, r, u);
      ursiz1 = REAL_RSQRT(DOT_(ur, ur));
      NORMAL_(ur, ursiz1);

      real us[3], ussiz1;
      CROSS_(us, s, u);
      ussiz1 = REAL_RSQRT(DOT_(us, us));
      NORMAL_(us, ussiz1);

      real vs[3], vssiz1;
      CROSS_(vs, s, v);
      vssiz1 = REAL_RSQRT(DOT_(vs, vs));
      NORMAL_(vs, vssiz1);

      real ws[3], wssiz1;
      CROSS_(ws, s, w);
      wssiz1 = REAL_RSQRT(DOT_(ws, ws));
      NORMAL_(ws, wssiz1);

      // get sine and cosine of angles between the rotation axes

      real urcos = DOT_(u, r);
      real ursin1 = REAL_RSQRT(1 - urcos * urcos);
      real vscos = DOT_(v, s);
      real vssin = REAL_SQRT(1 - vscos * vscos);
      real wscos = DOT_(w, s);
      real wssin = REAL_SQRT(1 - wscos * wscos);

      // compute the projection of v and w onto the ru-plane

      real t1[3], t1siz1;
      real t2[3], t2siz1;

      t1[0] = v[0] - s[0] * vscos;
      t1[1] = v[1] - s[1] * vscos;
      t1[2] = v[2] - s[2] * vscos;
      t1siz1 = REAL_RSQRT(DOT_(t1, t1));
      NORMAL_(t1, t1siz1);

      t2[0] = w[0] - s[0] * wscos;
      t2[1] = w[1] - s[1] * wscos;
      t2[2] = w[2] - s[2] * wscos;
      t2siz1 = REAL_RSQRT(DOT_(t2, t2));
      NORMAL_(t2, t2siz1);

      real ut1cos = DOT_(u, t1);
      real ut1sin = REAL_SQRT(1 - ut1cos * ut1cos);
      real ut2cos = DOT_(u, t2);
      real ut2sin = REAL_SQRT(1 - ut2cos * ut2cos);
      real _1_ut1sin_ut2sin = REAL_RECIP(ut1sin + ut2sin);

      // negative of dot product of torque with unit vectors gives result of
      // infinitesimal rotation along these vectors

      real dphidr = -DOT_(trq, r);
      real dphids = -DOT_(trq, s);

      for (int j = 0; j < 3; ++j) {
        real du = ur[j] * dphidr * usiz1 * ursin1 + us[j] * dphids * usiz1;
        real dv =
            (vssin * s[j] - vscos * t1[j]) * dphidu * vsiz1 * _1_ut1sin_ut2sin;
        real dw =
            (wssin * s[j] - wscos * t2[j]) * dphidu * wsiz1 * _1_ut1sin_ut2sin;
        #pragma acc atomic update
        de[j][ia] += du;
        #pragma acc atomic update
        de[j][ic] += dv;
        #pragma acc atomic update
        de[j][id] += dw;
        #pragma acc atomic update
        de[j][ib] -= (du + dv + dw);
        if_constexpr(DO_V) {
          frcz[j] += du;
          frcx[j] += dv;
          frcy[j] += dw;
        }
      }
    } else if (axetyp == pole_3_fold) {
      real p[3], psiz1;
      p[0] = u[0] + v[0] + w[0];
      p[1] = u[1] + v[1] + w[1];
      p[2] = u[2] + v[2] + w[2];
      psiz1 = REAL_RSQRT(DOT_(p, p));
      NORMAL_(p, psiz1);

      real wpcos = DOT_(w, p);
      real upcos = DOT_(u, p);
      real vpcos = DOT_(v, p);

      real r[3], rsiz1;
      real del[3], delsiz1;
      real eps[3], dphidr, dphiddel;

      ADD_(r, u, v);
      rsiz1 = REAL_RSQRT(DOT_(r, r));
      NORMAL_(r, rsiz1);
      real rwcos = DOT_(r, w);
      real rwsin1 = REAL_RSQRT(1 - rwcos * rwcos);
      dphidr = -DOT_(trq, r);
      CROSS_(del, r, w);
      delsiz1 = REAL_RSQRT(DOT_(del, del));
      NORMAL_(del, delsiz1);
      dphiddel = -DOT_(trq, del);
      CROSS_(eps, del, w);
      for (int j = 0; j < 3; ++j) {
        real dw = del[j] * dphidr * wsiz1 * rwsin1 +
            eps[j] * dphiddel * wpcos * wsiz1 * psiz1;
        #pragma acc atomic update
        de[j][id] += dw;
        #pragma acc atomic update
        de[j][ib] -= dw;
        if_constexpr(DO_V) frcy[j] += dw;
      }

      ADD_(r, v, w);
      rsiz1 = REAL_RSQRT(DOT_(r, r));
      NORMAL_(r, rsiz1);
      real rucos = DOT_(r, u);
      real rusin1 = REAL_RSQRT(1 - rucos * rucos);
      dphidr = -DOT_(trq, r);
      CROSS_(del, r, u);
      delsiz1 = REAL_RSQRT(DOT_(del, del));
      NORMAL_(del, delsiz1);
      dphiddel = -DOT_(trq, del);
      CROSS_(eps, del, u);
      for (int j = 0; j < 3; ++j) {
        real du = del[j] * dphidr * usiz1 * rusin1 +
            eps[j] * dphiddel * upcos * usiz1 * psiz1;
        #pragma acc atomic update
        de[j][ia] += du;
        #pragma acc atomic update
        de[j][ib] -= du;
        if_constexpr(DO_V) frcz[j] += du;
      }

      ADD_(r, u, w);
      rsiz1 = REAL_RSQRT(DOT_(r, r));
      NORMAL_(r, rsiz1);
      real rvcos = DOT_(r, v);
      real rvsin1 = REAL_RSQRT(1 - rvcos * rvcos);
      dphidr = -DOT_(trq, r);
      CROSS_(del, r, v);
      delsiz1 = REAL_RSQRT(DOT_(del, del));
      NORMAL_(del, delsiz1);
      dphiddel = -DOT_(trq, del);
      CROSS_(eps, del, v);
      for (int j = 0; j < 3; ++j) {
        real dv = del[j] * dphidr * vsiz1 * rvsin1 +
            eps[j] * dphiddel * vpcos * vsiz1 * psiz1;
        #pragma acc atomic update
        de[j][ic] += dv;
        #pragma acc atomic update
        de[j][ib] -= dv;
        if_constexpr(DO_V) frcx[j] += dv;
      }
    }

    if_constexpr(DO_V) {
      const int iaz = (ia == -1) ? i : ia;
      const int iax = (ic == -1) ? i : ic;
      const int iay = (id == -1) ? i : id;

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

      int offv = (i & (bufsize - 1)) * 16;
      atomic_add_value(gpu_vir, vxx, offv + 0);
      atomic_add_value(gpu_vir, vxy, offv + 1);
      atomic_add_value(gpu_vir, vxz, offv + 2);
      atomic_add_value(gpu_vir, vxy, offv + 3);
      atomic_add_value(gpu_vir, vyy, offv + 4);
      atomic_add_value(gpu_vir, vyz, offv + 5);
      atomic_add_value(gpu_vir, vxz, offv + 6);
      atomic_add_value(gpu_vir, vyz, offv + 7);
      atomic_add_value(gpu_vir, vzz, offv + 8);
    } // end if_constexpr(DO_V)
  }   // end for (int i)
}

void torque(int vers) {
  if (!use_elec())
    return;

  if (vers & calc::virial) {
    torque_tmpl<1>(vir_trq_handle);
  } else if (vers & calc::grad) {
    torque_tmpl<0>(-1);
  }
}
TINKER_NAMESPACE_END

#undef ADD_
#undef SUB_
#undef DOT_
#undef CROSS_
#undef NORMAL_
