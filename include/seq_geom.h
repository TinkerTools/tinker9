#pragma once
#include "add.h"
#include "mathfunc.h"
#include "seq_def.h"


namespace tinker {
#pragma acc routine seq
template <class Ver>
SEQ_CUDA
void dk_geom_group(
   real& restrict e, real& restrict vxx, real& restrict vyx, real& restrict vzx,
   real& restrict vyy, real& restrict vzy, real& restrict vzz,

   grad_prec* restrict degx, grad_prec* restrict degy, grad_prec* restrict degz,

   int i, const int (*restrict igfix)[2], const real (*restrict gfix)[3],

   const real* restrict x, const real* restrict y, const real* restrict z,
   const mass_prec* restrict mass, const int* restrict molec,
   const int (*restrict igrp)[2], const int* restrict kgrp,
   const mass_prec* restrict grpmass, TINKER_IMAGE_PARAMS)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   int ia = igfix[i][0];
   int ib = igfix[i][1];
   int ja1 = igrp[ia][0];
   int ja2 = igrp[ia][1];
   int jb1 = igrp[ib][0];
   int jb2 = igrp[ib][1];


   real xacm = 0;
   real yacm = 0;
   real zacm = 0;
   #pragma acc loop seq
   for (int j = ja1; j < ja2; ++j) {
      int k = kgrp[j];
      real weigh = mass[k];
      xacm += x[k] * weigh;
      yacm += y[k] * weigh;
      zacm += z[k] * weigh;
   }
   real weigha = REAL_MAX((real)1, grpmass[ia]);
   weigha = REAL_RECIP(weigha);


   real xbcm = 0;
   real ybcm = 0;
   real zbcm = 0;
   #pragma acc loop seq
   for (int j = jb1; j < jb2; ++j) {
      int k = kgrp[j];
      real weigh = mass[k];
      xbcm += x[k] * weigh;
      ybcm += y[k] * weigh;
      zbcm += z[k] * weigh;
   }
   real weighb = REAL_MAX((real)1, grpmass[ib]);
   weighb = REAL_RECIP(weighb);


   real xr = xacm * weigha - xbcm * weighb;
   real yr = yacm * weigha - ybcm * weighb;
   real zr = zacm * weigha - zbcm * weighb;


   bool intermol = molec[kgrp[ja1]] != molec[kgrp[jb1]];
   if (intermol)
      image(xr, yr, zr);


   real r = REAL_SQRT(xr * xr + yr * yr + zr * zr);
   real force = gfix[i][0];
   real gf1 = gfix[i][1];
   real gf2 = gfix[i][2];
   real target = (r < gf1 ? gf1 : (r > gf2 ? gf2 : r));
   real dt = r - target;


   if CONSTEXPR (do_e) {
      real dt2 = dt * dt;
      e = force * dt2;
   }
   if CONSTEXPR (do_g) {
      real rinv = (r == 0 ? 1 : REAL_RECIP(r));
      real de = 2 * force * dt * rinv;
      real dedx = de * xr;
      real dedy = de * yr;
      real dedz = de * zr;

      #pragma acc loop seq
      for (int j = ja1; j < ja2; ++j) {
         int k = kgrp[j];
         real ratio = mass[k] * weigha;
         atomic_add(dedx * ratio, degx, k);
         atomic_add(dedy * ratio, degy, k);
         atomic_add(dedz * ratio, degz, k);
      }
      #pragma acc loop seq
      for (int j = jb1; j < jb2; ++j) {
         int k = kgrp[j];
         real ratio = mass[k] * weighb;
         atomic_add(-dedx * ratio, degx, k);
         atomic_add(-dedy * ratio, degy, k);
         atomic_add(-dedz * ratio, degz, k);
      }
      if CONSTEXPR (do_v) {
         vxx = xr * dedx;
         vyx = yr * dedx;
         vzx = zr * dedx;
         vyy = yr * dedy;
         vzy = zr * dedy;
         vzz = zr * dedz;
      }
   }
}


#pragma acc routine seq
template <class Ver>
SEQ_CUDA
void dk_geom_distance(
   real& restrict e, real& restrict vxx, real& restrict vyx, real& restrict vzx,
   real& restrict vyy, real& restrict vzy, real& restrict vzz,

   grad_prec* restrict degx, grad_prec* restrict degy, grad_prec* restrict degz,

   int i, const int (*restrict idfix)[2], const real (*restrict dfix)[3],

   const real* restrict x, const real* restrict y, const real* restrict z,
   const int* restrict molec, TINKER_IMAGE_PARAMS)
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;


   int ia = idfix[i][0];
   int ib = idfix[i][1];
   real force = dfix[i][0];
   real df1 = dfix[i][1];
   real df2 = dfix[i][2];

   real xr = x[ia] - x[ib];
   real yr = y[ia] - y[ib];
   real zr = z[ia] - z[ib];
   bool intermol = molec[ia] != molec[ib];
   if (intermol)
      image(xr, yr, zr);
   real r = REAL_SQRT(xr * xr + yr * yr + zr * zr);
   real target = r;
   if (r < df1)
      target = df1;
   if (r > df2)
      target = df2;
   real dt = r - target;


   if CONSTEXPR (do_e) {
      real dt2 = dt * dt;
      e = force * dt2;
   }
   if CONSTEXPR (do_g) {
      if (r == 0)
         r = 1;
      real de = 2 * force * dt * REAL_RECIP(r);
      real dedx = de * xr;
      real dedy = de * yr;
      real dedz = de * zr;
      atomic_add(dedx, degx, ia);
      atomic_add(dedy, degy, ia);
      atomic_add(dedz, degz, ia);
      atomic_add(-dedx, degx, ib);
      atomic_add(-dedy, degy, ib);
      atomic_add(-dedz, degz, ib);
      if CONSTEXPR (do_v) {
         vxx = xr * dedx;
         vyx = yr * dedx;
         vzx = zr * dedx;
         vyy = yr * dedy;
         vzy = zr * dedy;
         vzz = zr * dedz;
      }
   }
}
}
