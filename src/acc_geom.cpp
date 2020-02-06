#include "add.h"
#include "box.h"
#include "e_geom.h"
#include "group.h"
#include "md.h"
#include "molecule.h"
#include "named_struct.h"
#include "seq_image.h"
#include <cassert>


TINKER_NAMESPACE_BEGIN
template <class Ver>
void egeom_acc1()
{
   constexpr int do_e = Ver::e;
   constexpr int do_g = Ver::g;
   constexpr int do_v = Ver::v;

   auto bufsize = buffer_size();

   const auto* molec = molecule.molecule;
   const auto* igrp = grp.igrp;
   const auto* kgrp = grp.kgrp;
   const auto* grpmass = grp.grpmass;

   #pragma acc kernels async deviceptr(x,y,z,gx,gy,gz,box,mass,molec,\
               igrp,kgrp,grpmass,igfix,gfix,\
               eg,vir_eg)
   {
      // group restraints
      #pragma acc loop independent
      for (int i = 0; i < ngfix; ++i) {
         int offset = i & (bufsize - 1);
         int ia = igfix[i][0];
         int ib = igfix[i][1];
         int ja1 = igrp[ia][0];
         int ja2 = igrp[ia][1];
         int jb1 = igrp[ib][0];
         int jb2 = igrp[ib][1];

         real xacm = 0;
         real yacm = 0;
         real zacm = 0;
         #pragma acc loop independent reduction(+:xacm,yacm,zacm)
         for (int j = ja1; j < ja2; ++j) {
            int k = kgrp[j];
            real weigh = mass[k];
            xacm += x[k] * weigh;
            yacm += y[k] * weigh;
            zacm += z[k] * weigh;
         }
         real weigha = REAL_MAX(1, grpmass[ia]);
         weigha = REAL_RECIP(weigha);

         real xbcm = 0;
         real ybcm = 0;
         real zbcm = 0;
         #pragma acc loop independent reduction(+:xbcm,ybcm,zbcm)
         for (int j = jb1; j < jb2; ++j) {
            int k = kgrp[j];
            real weigh = mass[k];
            xbcm += x[k] * weigh;
            ybcm += y[k] * weigh;
            zbcm += z[k] * weigh;
         }
         real weighb = REAL_MAX(1, grpmass[ib]);
         weighb = REAL_RECIP(weighb);

         real xr = xacm * weigha - xbcm * weighb;
         real yr = yacm * weigha - ybcm * weighb;
         real zr = zacm * weigha - zbcm * weighb;

         bool intermol = molec[kgrp[ja1]] != molec[kgrp[jb1]];
         if (intermol)
            image(xr, yr, zr, box);

         real r = REAL_SQRT(xr * xr + yr * yr + zr * zr);
         real force = gfix[i][0];
         real gf1 = gfix[i][1];
         real gf2 = gfix[i][2];
         real target = (r < gf1 ? gf1 : (r > gf2 ? gf2 : r));
         real dt = r - target;

         if CONSTEXPR (do_e) {
            real dt2 = dt * dt;
            real e = force * dt2;
            atomic_add(e, eg, offset);
         }
         if CONSTEXPR (do_g) {
            real rinv = (r == 0 ? 1 : REAL_RECIP(r));
            real de = 2 * force * dt * rinv;
            real dedx = de * xr;
            real dedy = de * yr;
            real dedz = de * zr;

            #pragma acc loop independent
            for (int j = ja1; j < ja2; ++j) {
               int k = kgrp[j];
               real ratio = mass[k] * weigha;
               atomic_add(dedx * ratio, gx, k);
               atomic_add(dedy * ratio, gy, k);
               atomic_add(dedz * ratio, gz, k);
            }

            #pragma acc loop independent
            for (int j = jb1; j < jb2; ++j) {
               int k = kgrp[j];
               real ratio = mass[k] * weighb;
               atomic_add(-dedx * ratio, gx, k);
               atomic_add(-dedy * ratio, gy, k);
               atomic_add(-dedz * ratio, gz, k);
            }

            if CONSTEXPR (do_v) {
               real vxx = xr * dedx;
               real vyx = yr * dedx;
               real vzx = zr * dedx;
               real vyy = yr * dedy;
               real vzy = zr * dedy;
               real vzz = zr * dedz;

               atomic_add(vxx, vyx, vzx, vyy, vzy, vzz, vir_eg, offset);
            } // end if (do_v)
         }    // end if (do_g)
      }       // end ngfix loop
   }
}


void egeom_acc(int vers)
{
   if (vers == calc::v0 || vers == calc::v3)
      egeom_acc1<EnergyVersion0>();
   else if (vers == calc::v1)
      egeom_acc1<EnergyVersion1>();
   else if (vers == calc::v4)
      egeom_acc1<EnergyVersion4>();
   else if (vers == calc::v5)
      egeom_acc1<EnergyVersion5>();
   else if (vers == calc::v6)
      egeom_acc1<EnergyVersion6>();
   else
      assert(false);
}
TINKER_NAMESPACE_END
