#include "acc_seq.h"
#include "couple.h"
#include "gpu/e_vdw.h"
#include "md.h"
#include "nblist.h"
#include <vector>

// TODO: test lj, buck, mm3hb, gauss, and mutant
// TODO: add vdw correction

TINKER_NAMESPACE_BEGIN
void evdw_reduce_xyz() {
  #pragma acc parallel loop independent\
              deviceptr(x,y,z,ired,kred,xred,yred,zred)
  for (int i = 0; i < n; ++i) {
    int iv = ired[i];
    real rdn = kred[i];
    xred[i] = rdn * (x[i] - x[iv]) + x[iv];
    yred[i] = rdn * (y[i] - y[iv]) + y[iv];
    zred[i] = rdn * (z[i] - z[iv]) + z[iv];
  }
}

template <int USE, evdw_t VDWTYP>
void evdw_tmpl() {
  constexpr int do_e = USE & calc::energy;
  constexpr int do_a = USE & calc::analyz;
  constexpr int do_g = USE & calc::grad;
  constexpr int do_v = USE & calc::virial;
  static_assert(do_v ? do_g : true, "");
  static_assert(do_a ? do_e : true, "");

  const real cut = vdw_switch_cut;
  const real off = vdw_switch_off;
  const real cut2 = cut * cut;
  const real off2 = off * off;
  const int maxnlst = vlist_unit.obj().maxnlst;
  const NBList* vlst = vlist_unit.deviceptr();

  const auto* coupl = couple_unit.deviceptr();

  static std::vector<real> vscalebuf;
  vscalebuf.resize(n, 1);
  real* vscale = vscalebuf.data();

  #pragma acc parallel loop independent\
              deviceptr(x,y,z,gx,gy,gz,box,coupl,vlst,\
                        ired,kred,xred,yred,zred,\
                        jvdw,njvdw,radmin,epsilon,vlam,\
                        ev,nev,vir_ev)\
              firstprivate(vscale[0:n])
  for (int i = 0; i < n; ++i) {
    const int n12i = coupl->n12[i];
    const int n13i = coupl->n13[i];
    const int n14i = coupl->n14[i];
    const int n15i = coupl->n15[i];
    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j)
      vscale[coupl->i12[i][j]] = v2scale;
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j)
      vscale[coupl->i13[i][j]] = v3scale;
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j)
      vscale[coupl->i14[i][j]] = v4scale;
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j)
      vscale[coupl->i15[i][j]] = v5scale;

    int iv = ired[i];
    real redi = kred[i];
    real rediv = 1 - redi;
    int it = jvdw[i];
    real xi, yi, zi;
    xi = xred[i];
    yi = yred[i];
    zi = zred[i];
    real lam1 = vlam[i];

    int base_it = it * (*njvdw);
    int nvlsti = vlst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop independent
    for (int kk = 0; kk < nvlsti; ++kk) {
      int k = vlst->lst[base + kk];
      int kv = ired[k];
      int kt = jvdw[k];
      real xr, yr, zr;
      xr = xi - xred[k];
      yr = yi - yred[k];
      zr = zi - zred[k];
      real vlambda = vlam[k];

      if (vcouple == vcouple_decouple) {
        vlambda = (lam1 == vlambda ? 1 : (lam1 < vlambda ? lam1 : vlambda));
      } else if (vcouple == vcouple_annihilate) {
        vlambda = (lam1 < vlambda ? lam1 : vlambda);
      }

      image(xr, yr, zr, box);
      real rik2 = xr * xr + yr * yr + zr * zr;
      if (rik2 <= off2) {
        real rik = REAL_SQRT(rik2);
        real rv = radmin[base_it + kt];
        real eps = epsilon[base_it + kt] * vscale[k];

        real e;
        real de;
        if_constexpr(VDWTYP & vdw_hal) {
          real rho = rik * REAL_RECIP(rv);
          real rho6 = REAL_POW(rho, 6);
          real rho7 = rho6 * rho;
          eps *= REAL_POW(vlambda, scexp);
          real scal = scalpha * REAL_POW(1 - vlambda, 2);
          real s1 = REAL_RECIP(scal + REAL_POW(rho + dhal, 7));
          real s2 = REAL_RECIP(scal + rho7 + ghal);
          real t1 = REAL_POW(1 + dhal, 7) * s1;
          real t2 = (1 + ghal) * s2;
          e = eps * t1 * (t2 - 2);
          if_constexpr(do_g) {
            real dt1drho = -7 * REAL_POW(rho + dhal, 6) * t1 * s1;
            real dt2drho = -7 * rho6 * t2 * s2;
            de = eps * (dt1drho * (t2 - 2) + t1 * dt2drho) * REAL_RECIP(rv);
          }
        }

        if (rik2 > cut2) {
          real taper, dtaper;
          switch_taper5<do_g>(rik, cut, off, taper, dtaper);
          if_constexpr(do_g) de = e * dtaper + de * taper;
          if_constexpr(do_e) e = e * taper;
        }

        // Increment the energy, gradient, and virial.

        if_constexpr(do_e) {
          #pragma acc atomic update
          *ev += e;
          if_constexpr(do_a) {
            if (e != 0) {
              #pragma acc atomic update
              *nev += 1;
            }
          }
        }

        real dedx, dedy, dedz;
        if_constexpr(do_g) {
          de *= REAL_RECIP(rik);
          dedx = de * xr;
          dedy = de * yr;
          dedz = de * zr;

          #pragma acc atomic update
          gx[i] += dedx * redi;
          #pragma acc atomic update
          gy[i] += dedy * redi;
          #pragma acc atomic update
          gz[i] += dedz * redi;
          #pragma acc atomic update
          gx[iv] += dedx * rediv;
          #pragma acc atomic update
          gy[iv] += dedy * rediv;
          #pragma acc atomic update
          gz[iv] += dedz * rediv;

          real redk = kred[k];
          real redkv = 1 - redk;
          #pragma acc atomic update
          gx[k] -= dedx * redk;
          #pragma acc atomic update
          gy[k] -= dedy * redk;
          #pragma acc atomic update
          gz[k] -= dedz * redk;
          #pragma acc atomic update
          gx[kv] -= dedx * redkv;
          #pragma acc atomic update
          gy[kv] -= dedy * redkv;
          #pragma acc atomic update
          gz[kv] -= dedz * redkv;
        }

        if_constexpr(do_v) {
          real vxx = xr * dedx;
          real vyx = yr * dedx;
          real vzx = zr * dedx;
          real vyy = yr * dedy;
          real vzy = zr * dedy;
          real vzz = zr * dedz;

          #pragma acc atomic update
          vir_ev[_xx] += vxx;
          #pragma acc atomic update
          vir_ev[_yx] += vyx;
          #pragma acc atomic update
          vir_ev[_zx] += vzx;
          #pragma acc atomic update
          vir_ev[_xy] += vyx;
          #pragma acc atomic update
          vir_ev[_yy] += vyy;
          #pragma acc atomic update
          vir_ev[_zy] += vzy;
          #pragma acc atomic update
          vir_ev[_xz] += vzx;
          #pragma acc atomic update
          vir_ev[_yz] += vzy;
          #pragma acc atomic update
          vir_ev[_zz] += vzz;
        }
      }
    } // end for (int kk)

    #pragma acc loop independent
    for (int j = 0; j < n12i; ++j)
      vscale[coupl->i12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n13i; ++j)
      vscale[coupl->i13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n14i; ++j)
      vscale[coupl->i14[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < n15i; ++j)
      vscale[coupl->i15[i][j]] = 1;
  } // end for (int i)
}

#define TINKER_EVDW_IMPL_(typ)                                                 \
  void evdw_##typ##_acc_impl_(int vers) {                                      \
    evdw_reduce_xyz();                                                         \
    if (vers == calc::v0)                                                      \
      evdw_tmpl<calc::v0, vdw_##typ>();                                        \
    else if (vers == calc::v1)                                                 \
      evdw_tmpl<calc::v1, vdw_##typ>();                                        \
    else if (vers == calc::v3)                                                 \
      evdw_tmpl<calc::v3, vdw_##typ>();                                        \
    else if (vers == calc::v4)                                                 \
      evdw_tmpl<calc::v4, vdw_##typ>();                                        \
    else if (vers == calc::v5)                                                 \
      evdw_tmpl<calc::v5, vdw_##typ>();                                        \
    else if (vers == calc::v6)                                                 \
      evdw_tmpl<calc::v6, vdw_##typ>();                                        \
  }
TINKER_EVDW_IMPL_(lj);
TINKER_EVDW_IMPL_(buck);
TINKER_EVDW_IMPL_(mm3hb);
TINKER_EVDW_IMPL_(hal);
TINKER_EVDW_IMPL_(gauss);
#undef TINKER_EVDW_IMPL_

TINKER_NAMESPACE_END
