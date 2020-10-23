#include "add.h"
#include "cflux.h"
#include "couple.h"
#include "eangle.h"
#include "ebond.h"
#include "elec.h"
#include "mathfunc.h"
#include "md.h"
#include "seq_bsplgen.h"
#include "tool/energy_buffer.h"
#include "tool/gpu_card.h"

namespace tinker {
#pragma acc routine seq
real dot_vect(const real* restrict a, const real* restrict b)
{
   return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

void bnd_dcflux()
{
   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,ibnd,bflx,pot,\
               decfx,decfy,decfz)
   for (int i = 0; i < nbond; ++i) {
      int ia = ibnd[i][0];
      int ib = ibnd[i][1];
      real pb = bflx[i];
      real xab = x[ia] - x[ib];
      real yab = y[ia] - y[ib];
      real zab = z[ia] - z[ib];

      real rab = REAL_SQRT(xab * xab + yab * yab + zab * zab);
      real dpot = pot[ia] - pot[ib];
      pb = pb / rab;

      real fx = dpot * pb * xab;
      real fy = dpot * pb * yab;
      real fz = dpot * pb * zab;

      atomic_add(-fx, decfx, ia);
      atomic_add(-fy, decfy, ia);
      atomic_add(-fz, decfz, ia);
      atomic_add(fx, decfx, ib);
      atomic_add(fy, decfy, ib);
      atomic_add(fz, decfz, ib);
   }
}

void ang_dcflux()
{
   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,iang,aflx,abflx,pot,\
               decfx,decfy,decfz)
   for (int i = 0; i < nangle; ++i) {
      int ia = iang[i][0];
      int ib = iang[i][1];
      int ic = iang[i][2];
      real pa1 = aflx[i][0];
      real pa2 = aflx[i][1];
      real pb1 = abflx[i][0];
      real pb2 = abflx[i][1];
      real xab = x[ia] - x[ib];
      real yab = y[ia] - y[ib];
      real zab = z[ia] - z[ib];
      real xcb = x[ic] - x[ib];
      real ycb = y[ic] - y[ib];
      real zcb = z[ic] - z[ib];

      real rba2 = xab * xab + yab * yab + zab * zab;
      real rba = REAL_SQRT(rba2);
      real rba3 = rba2 * rba;
      real dba[3];
      dba[0] = xab;
      dba[1] = yab;
      dba[2] = zab;

      real rbc2 = xcb * xcb + ycb * ycb + zcb * zcb;
      real rbc = REAL_SQRT(rbc2);
      real rbc3 = rbc2 * rbc;
      real dbc[3];
      dbc[0] = xcb;
      dbc[1] = ycb;
      dbc[2] = zcb;

      real dpota = pot[ia] - pot[ib];
      real dpotc = pot[ic] - pot[ib];
      pb1 *= dpota;
      pb2 *= dpotc;

      real c1 = pb2 / rba;
      real c2 = pb1 / rbc;

      real fax = c1 * xab;
      real fay = c1 * yab;
      real faz = c1 * zab;
      real fcx = c2 * xcb;
      real fcy = c2 * ycb;
      real fcz = c2 * zcb;
      real fbx = -(fax + fcx);
      real fby = -(fay + fcy);
      real fbz = -(faz + fcz);

      real dot = dot_vect(dba, dbc);
      real term = -rba * rbc / REAL_SQRT(rba2 * rbc2 - dot * dot);
      real fterm = term * (dpota * pa1 + dpotc * pa2);

      c1 = 1 / (rba * rbc);
      c2 = dot / (rba3 * rbc);
      real c3 = dot / (rbc3 * rba);

      real fax2 = fterm * (c1 * xcb - c2 * xab);
      real fay2 = fterm * (c1 * ycb - c2 * yab);
      real faz2 = fterm * (c1 * zcb - c2 * zab);
      real fcx2 = fterm * (c1 * xab - c3 * xcb);
      real fcy2 = fterm * (c1 * yab - c3 * ycb);
      real fcz2 = fterm * (c1 * zab - c3 * zcb);
      real fbx2 = -(fax2 + fcx2);
      real fby2 = -(fay2 + fcy2);
      real fbz2 = -(faz2 + fcz2);

      atomic_add((fax + fax2), decfx, ia);
      atomic_add((fay + fay2), decfy, ia);
      atomic_add((faz + faz2), decfz, ia);
      atomic_add((fbx + fbx2), decfx, ib);
      atomic_add((fby + fby2), decfy, ib);
      atomic_add((fbz + fbz2), decfz, ib);
      atomic_add((fcx + fcx2), decfx, ic);
      atomic_add((fcy + fcy2), decfy, ic);
      atomic_add((fcz + fcz2), decfz, ic);
   }
}

template <int DO_V>
void dcflux_acc1(grad_prec* restrict gx, grad_prec* restrict gy,
                 grad_prec* restrict gz, virial_buffer restrict vir)
{
   auto bufsize = buffer_size();

   #pragma acc parallel loop independent async\
               deviceptr(decfx,decfy,decfz)
   for (int i = 0; i < n; ++i) {
      decfx[i] = 0;
      decfy[i] = 0;
      decfz[i] = 0;
   }

   bnd_dcflux();
   ang_dcflux();

   #pragma acc parallel loop independent async\
               deviceptr(x,y,z,decfx,decfy,decfz,\
               gx,gy,gz,vir)
   for (int i = 0; i < n; ++i) {
      atomic_add(decfx[i], gx, i);
      atomic_add(decfy[i], gy, i);
      atomic_add(decfz[i], gz, i);

      if CONSTEXPR (DO_V) {
         real vxx = x[i] * decfx[i];
         real vyx = y[i] * decfx[i];
         real vzx = z[i] * decfx[i];
         real vyy = y[i] * decfy[i];
         real vzy = z[i] * decfy[i];
         real vzz = z[i] * decfz[i];
         atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir, i & (bufsize - 1));
      }
   }
}

void dcflux_acc(int vers, grad_prec* gx, grad_prec* gy, grad_prec* gz,
                virial_buffer vir)
{
   if (vers & calc::virial)
      dcflux_acc1<1>(gx, gy, gz, vir);
   else
      dcflux_acc1<0>(gx, gy, gz, nullptr);
}
}
