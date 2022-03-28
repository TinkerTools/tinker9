#include "add.h"
#include "ff/elec.h"
#include "math/inc.h"
#include "md/inc.h"
#include "mod/elecamoeba.h"

namespace tinker {
namespace {
#pragma acc routine seq
real torque_dot(const real* restrict a, const real* restrict b)
{
   return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

#pragma acc routine seq
void torque_cross(real* restrict ans, const real* restrict u, const real* restrict v)
{
   ans[0] = u[1] * v[2] - u[2] * v[1];
   ans[1] = u[2] * v[0] - u[0] * v[2];
   ans[2] = u[0] * v[1] - u[1] * v[0];
}

#pragma acc routine seq
void torque_normal(real* restrict a, real _1_na)
{
   a[0] *= _1_na;
   a[1] *= _1_na;
   a[2] *= _1_na;
}
}

/**
 * This commit `09beaa8b0215913c1b658b0b2d8ecdf97bda8c15` had a weird problem
 * in the OpenACC torque kernel. It had a strange segfault when calling torque
 * in test "NaCl-2" after "Local-Frame-3" having been called. I think the
 * problem might be on PGI's side, especially when I don't see any problems
 * reported by `cuda-memcheck`. Run this command will see the segfault
 * \code
 * GPU_PACKAGE=OPENACC PGI_ACC_NOTIFY=31 PGI_ACC_DEBUG=1 \
 * ./all.tests local-frame-3,nacl-2 \
 * -d yes -a
 * \endcode
 *
 * This kernel has a lot of float[3] arrays private to each thread, and the PGI
 * compiler decided not to use the "local memory" for two of them, ws[3] and
 * vs[3], which resulted in 24 bytes memory copyin at the beginning and copyout
 * at the end of the routine. When I have rewritten the torque kernel to force
 * everything in the "local memory", the segfault would disappear.
 */
template <int DO_V>
void torque_acc1(VirialBuffer gpu_vir, grad_prec* gx, grad_prec* gy, grad_prec* gz)
{
   auto bufsize = bufferSize();

   real trq[3], frcz[3], frcx[3], frcy[3];
   real deia[3], deib[3], deic[3], deid[3];
   real u[3], v[3], w[3];
   // uv uw: z only, z then x
   // vw:    z only, z then x, bisector
   real uv[3], uw[3], vw[3];
   // z bisector
   real r[3], s[3], ur[3], us[3], vs[3], ws[3], t1[3], t2[3];
   // 3 fold
   real p[3], del[3], eps[3];
   #pragma acc parallel loop independent async private(\
               trq[3],frcz[3],frcx[3],frcy[3],\
               deia[3],deib[3],deic[3],deid[3],\
               u[3],v[3],w[3],uv[3],uw[3],vw[3],\
               r[3],s[3],ur[3],us[3],vs[3],ws[3],t1[3],t2[3],\
               p[3],del[3],eps[3])\
               deviceptr(x,y,z,gx,gy,gz,zaxis,trqx,trqy,trqz,gpu_vir)
   for (int i = 0; i < n; ++i) {
      const int axetyp = zaxis[i].polaxe;
      if (axetyp == pole_none)
         continue;

      // get the local frame type and the frame-defining atoms
      const int ia = zaxis[i].zaxis;
      const int ib = i;
      const int ic = zaxis[i].xaxis;
      const int id = INT_ABS(zaxis[i].yaxis) - 1;

      // construct the three rotation axes for the local frame
      // u
      real usiz1;
      u[0] = x[ia] - x[ib];
      u[1] = y[ia] - y[ib];
      u[2] = z[ia] - z[ib];
      usiz1 = REAL_RSQRT(torque_dot(u, u));
      torque_normal(u, usiz1);

      // v
      real vsiz1;
      if (axetyp != pole_z_only) {
         v[0] = x[ic] - x[ib];
         v[1] = y[ic] - y[ib];
         v[2] = z[ic] - z[ib];
         vsiz1 = REAL_RSQRT(torque_dot(v, v));
      } else {
         int foo = !(REAL_ABS(u[0]) > (real)0.866);
         v[0] = (foo ? 1 : 0);
         v[1] = (foo ? 0 : 1);
         v[2] = 0;
         vsiz1 = 1;
      }
      torque_normal(v, vsiz1);

      // w = u cross v
      real wsiz1;
      if (axetyp == pole_z_bisect || axetyp == pole_3_fold) {
         w[0] = x[id] - x[ib];
         w[1] = y[id] - y[ib];
         w[2] = z[id] - z[ib];
      } else {
         torque_cross(w, u, v);
      }
      wsiz1 = REAL_RSQRT(torque_dot(w, w));
      torque_normal(w, wsiz1);

      // negative of dot product of torque with unit vectors gives result of
      // infinitesimal rotation along these vectors
      trq[0] = trqx[i];
      trq[1] = trqy[i];
      trq[2] = trqz[i];
      real dphidu = -torque_dot(trq, u);
      real dphidv = -torque_dot(trq, v);
      real dphidw = -torque_dot(trq, w);

      // force distribution
      deia[0] = 0;
      deia[1] = 0;
      deia[2] = 0;
      deib[0] = 0;
      deib[1] = 0;
      deib[2] = 0;
      deic[0] = 0;
      deic[1] = 0;
      deic[2] = 0;
      deid[0] = 0;
      deid[1] = 0;
      deid[2] = 0;

      // zero out force components on local frame-defining atoms
      frcz[0] = 0;
      frcz[1] = 0;
      frcz[2] = 0;
      frcx[0] = 0;
      frcx[1] = 0;
      frcx[2] = 0;
      frcy[0] = 0;
      frcy[1] = 0;
      frcy[2] = 0;

      if (axetyp == pole_z_only || axetyp == pole_z_then_x || axetyp == pole_bisector) {
         // find the perpendicular and angle for each pair of axes
         // v cross u
         real uvsiz1;
         torque_cross(uv, v, u);
         uvsiz1 = REAL_RSQRT(torque_dot(uv, uv));
         torque_normal(uv, uvsiz1);

         // w cross u
         real uwsiz1;
         torque_cross(uw, w, u);
         uwsiz1 = REAL_RSQRT(torque_dot(uw, uw));
         torque_normal(uw, uwsiz1);

         // get sine and cosine of angles between the rotation axes
         real uvcos = torque_dot(u, v);
         real uvsin1 = REAL_RSQRT(1 - uvcos * uvcos);

         if (axetyp == pole_z_only) {
            #pragma acc loop seq
            for (int j = 0; j < 3; ++j) {
               real du = uv[j] * dphidv * usiz1 * uvsin1 + uw[j] * dphidw * usiz1;
               deia[j] += du;
               deib[j] -= du;
               if CONSTEXPR (DO_V) {
                  frcz[j] += du;
               }
            }
         } else if (axetyp == pole_z_then_x) {
            #pragma acc loop seq
            for (int j = 0; j < 3; ++j) {
               real du = uv[j] * dphidv * usiz1 * uvsin1 + uw[j] * dphidw * usiz1;
               real dv = -uv[j] * dphidu * vsiz1 * uvsin1;
               deia[j] += du;
               deic[j] += dv;
               deib[j] -= (du + dv);
               if CONSTEXPR (DO_V) {
                  frcz[j] += du;
                  frcx[j] += dv;
               }
            }
         } else /* if (axetyp == pole_bisector) */ {
            // w cross v
            real vwsiz1;
            torque_cross(vw, w, v);
            vwsiz1 = REAL_RSQRT(torque_dot(vw, vw));
            torque_normal(vw, vwsiz1);

            #pragma acc loop seq
            for (int j = 0; j < 3; ++j) {
               real du = uv[j] * dphidv * usiz1 * uvsin1 + 0.5f * uw[j] * dphidw * usiz1;
               real dv = -uv[j] * dphidu * vsiz1 * uvsin1 + 0.5f * vw[j] * dphidw * vsiz1;
               deia[j] += du;
               deic[j] += dv;
               deib[j] -= (du + dv);
               if CONSTEXPR (DO_V) {
                  frcz[j] += du;
                  frcx[j] += dv;
               }
            }
         }
      } else if (axetyp == pole_z_bisect) {
         // build some additional axes needed for the Z-Bisect method
         real rsiz1;
         r[0] = v[0] + w[0];
         r[1] = v[1] + w[1];
         r[2] = v[2] + w[2];
         rsiz1 = REAL_RSQRT(torque_dot(r, r));

         real ssiz1;
         torque_cross(s, u, r);
         ssiz1 = REAL_RSQRT(torque_dot(s, s));

         torque_normal(r, rsiz1);
         torque_normal(s, ssiz1);

         // find the perpendicular and angle for each pair of axes
         real ursiz1;
         torque_cross(ur, r, u);
         ursiz1 = REAL_RSQRT(torque_dot(ur, ur));
         torque_normal(ur, ursiz1);

         real ussiz1;
         torque_cross(us, s, u);
         ussiz1 = REAL_RSQRT(torque_dot(us, us));
         torque_normal(us, ussiz1);

         real vssiz1;
         torque_cross(vs, s, v);
         vssiz1 = REAL_RSQRT(torque_dot(vs, vs));
         torque_normal(vs, vssiz1);

         real wssiz1;
         torque_cross(ws, s, w);
         wssiz1 = REAL_RSQRT(torque_dot(ws, ws));
         torque_normal(ws, wssiz1);

         // get sine and cosine of angles between the rotation axes
         real urcos = torque_dot(u, r);
         real ursin1 = REAL_RSQRT(1 - urcos * urcos);
         real vscos = torque_dot(v, s);
         real vssin = REAL_SQRT(1 - vscos * vscos);
         real wscos = torque_dot(w, s);
         real wssin = REAL_SQRT(1 - wscos * wscos);

         // compute the projection of v and w onto the ru-plane
         real t1siz1;
         t1[0] = v[0] - s[0] * vscos;
         t1[1] = v[1] - s[1] * vscos;
         t1[2] = v[2] - s[2] * vscos;
         t1siz1 = REAL_RSQRT(torque_dot(t1, t1));
         torque_normal(t1, t1siz1);

         real t2siz1;
         t2[0] = w[0] - s[0] * wscos;
         t2[1] = w[1] - s[1] * wscos;
         t2[2] = w[2] - s[2] * wscos;
         t2siz1 = REAL_RSQRT(torque_dot(t2, t2));
         torque_normal(t2, t2siz1);

         real ut1cos = torque_dot(u, t1);
         real ut1sin = REAL_SQRT(1 - ut1cos * ut1cos);
         real ut2cos = torque_dot(u, t2);
         real ut2sin = REAL_SQRT(1 - ut2cos * ut2cos);
         real _1_ut1sin_ut2sin = REAL_RECIP(ut1sin + ut2sin);

         // negative of dot product of torque with unit vectors gives result of
         // infinitesimal rotation along these vectors
         real dphidr = -torque_dot(trq, r);
         real dphids = -torque_dot(trq, s);

         #pragma acc loop seq
         for (int j = 0; j < 3; ++j) {
            real du = ur[j] * dphidr * usiz1 * ursin1 + us[j] * dphids * usiz1;
            real dv = (vssin * s[j] - vscos * t1[j]) * dphidu * vsiz1 * _1_ut1sin_ut2sin;
            real dw = (wssin * s[j] - wscos * t2[j]) * dphidu * wsiz1 * _1_ut1sin_ut2sin;
            deia[j] += du;
            deic[j] += dv;
            deid[j] += dw;
            deib[j] -= (du + dv + dw);
            if CONSTEXPR (DO_V) {
               frcz[j] += du;
               frcx[j] += dv;
               frcy[j] += dw;
            }
         }
      } else if (axetyp == pole_3_fold) {
         real psiz1;
         p[0] = u[0] + v[0] + w[0];
         p[1] = u[1] + v[1] + w[1];
         p[2] = u[2] + v[2] + w[2];
         psiz1 = REAL_RSQRT(torque_dot(p, p));
         torque_normal(p, psiz1);

         real wpcos = torque_dot(w, p);
         real upcos = torque_dot(u, p);
         real vpcos = torque_dot(v, p);

         real rsiz1;
         real delsiz1;
         real dphidr, dphiddel;

         r[0] = u[0] + v[0];
         r[1] = u[1] + v[1];
         r[2] = u[2] + v[2];
         rsiz1 = REAL_RSQRT(torque_dot(r, r));
         torque_normal(r, rsiz1);
         real rwcos = torque_dot(r, w);
         real rwsin1 = REAL_RSQRT(1 - rwcos * rwcos);
         dphidr = -torque_dot(trq, r);
         torque_cross(del, r, w);
         delsiz1 = REAL_RSQRT(torque_dot(del, del));
         torque_normal(del, delsiz1);
         dphiddel = -torque_dot(trq, del);
         torque_cross(eps, del, w);
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j) {
            real dw = del[j] * dphidr * wsiz1 * rwsin1 + eps[j] * dphiddel * wpcos * wsiz1 * psiz1;
            deid[j] += dw;
            deib[j] -= dw;
            if CONSTEXPR (DO_V) {
               frcy[j] += dw;
            }
         }

         r[0] = v[0] + w[0];
         r[1] = v[1] + w[1];
         r[2] = v[2] + w[2];
         rsiz1 = REAL_RSQRT(torque_dot(r, r));
         torque_normal(r, rsiz1);
         real rucos = torque_dot(r, u);
         real rusin1 = REAL_RSQRT(1 - rucos * rucos);
         dphidr = -torque_dot(trq, r);
         torque_cross(del, r, u);
         delsiz1 = REAL_RSQRT(torque_dot(del, del));
         torque_normal(del, delsiz1);
         dphiddel = -torque_dot(trq, del);
         torque_cross(eps, del, u);
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j) {
            real du = del[j] * dphidr * usiz1 * rusin1 + eps[j] * dphiddel * upcos * usiz1 * psiz1;
            deia[j] += du;
            deib[j] -= du;
            if CONSTEXPR (DO_V) {
               frcz[j] += du;
            }
         }

         r[0] = u[0] + w[0];
         r[1] = u[1] + w[1];
         r[2] = u[2] + w[2];
         rsiz1 = REAL_RSQRT(torque_dot(r, r));
         torque_normal(r, rsiz1);
         real rvcos = torque_dot(r, v);
         real rvsin1 = REAL_RSQRT(1 - rvcos * rvcos);
         dphidr = -torque_dot(trq, r);
         torque_cross(del, r, v);
         delsiz1 = REAL_RSQRT(torque_dot(del, del));
         torque_normal(del, delsiz1);
         dphiddel = -torque_dot(trq, del);
         torque_cross(eps, del, v);
         #pragma acc loop seq
         for (int j = 0; j < 3; ++j) {
            real dv = del[j] * dphidr * vsiz1 * rvsin1 + eps[j] * dphiddel * vpcos * vsiz1 * psiz1;
            deic[j] += dv;
            deib[j] -= dv;
            if CONSTEXPR (DO_V) {
               frcx[j] += dv;
            }
         }
      }

      atomic_add(deia[0], gx, ia);
      atomic_add(deia[1], gy, ia);
      atomic_add(deia[2], gz, ia);
      atomic_add(deib[0], gx, ib);
      atomic_add(deib[1], gy, ib);
      atomic_add(deib[2], gz, ib);
      if (axetyp != pole_z_only) {
         atomic_add(deic[0], gx, ic);
         atomic_add(deic[1], gy, ic);
         atomic_add(deic[2], gz, ic);
      }
      if (axetyp == pole_z_bisect || axetyp == pole_3_fold) {
         atomic_add(deid[0], gx, id);
         atomic_add(deid[1], gy, id);
         atomic_add(deid[2], gz, id);
      }

      if CONSTEXPR (DO_V) {
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
            (yix * frcx[0] + yiy * frcy[0] + yiz * frcz[0] + xix * frcx[1] + xiy * frcy[1] +
               xiz * frcz[1]);
         real vxz = 0.5f *
            (zix * frcx[0] + ziy * frcy[0] + ziz * frcz[0] + xix * frcx[2] + xiy * frcy[2] +
               xiz * frcz[2]);
         real vyy = yix * frcx[1] + yiy * frcy[1] + yiz * frcz[1];
         real vyz = 0.5f *
            (zix * frcx[1] + ziy * frcy[1] + ziz * frcz[1] + yix * frcx[2] + yiy * frcy[2] +
               yiz * frcz[2]);
         real vzz = zix * frcx[2] + ziy * frcy[2] + ziz * frcz[2];

         atomic_add(vxx, vxy, vxz, vyy, vyz, vzz, gpu_vir, i & (bufsize - 1));
      } // end if CONSTEXPR (DO_V)
   }    // end for (int i)
}

void torque_acc(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz)
{
   if (vers & calc::virial)
      torque_acc1<1>(vir_trq, gx, gy, gz);
   else if (vers & calc::grad)
      torque_acc1<0>(nullptr, gx, gy, gz);
}
}
