#include "group.h"
#include "md.h"
#include "tool/darray.h"
#include <tinker/detail/group.hh>
#include <tinker/detail/sizes.hh>

namespace tinker {
void group_data(rc_op op)
{
   if (op & rc_dealloc) {
      auto& st = grp;
      darray::deallocate(st.kgrp, st.grplist, st.igrp, st.grpmass, st.wgrp);
   }

   if (op & rc_alloc) {
      auto& st = grp;
      st.ngrp = group::ngrp;
      // see cluster.f
      // allocate (kgrp(n))
      // allocate (grplist(n))
      // allocate (igrp(2,0:maxgrp))
      // allocate (grpmass(0:maxgrp))
      // allocate (wgrp(0:maxgrp,0:maxgrp))
      darray::allocate(n, &st.kgrp, &st.grplist);
      darray::allocate(1 + st.ngrp, &st.igrp, &st.grpmass);
      darray::allocate((1 + st.ngrp) * (1 + st.ngrp), &st.wgrp);
   }

   if (op & rc_init) {
      // at most (n + 1) groups, e.g. NaCl ion pair
      // group 0: (NULL), igrp (fortran) 1,0, igrp (here) 0,0
      // group 1: Na+,    igrp (fortran) 1,1, igrp (here) 0,1
      // group 2: Cl-,    igrp (fortran) 2,2, igrp (here) 1,2

      auto& st = grp;
      std::vector<int> buf(2 * n + 2);

      for (int i = 0; i < n; ++i) {
         buf[i] = group::kgrp[i] - 1;
      }
      darray::copyin(g::q0, n, st.kgrp, buf.data());
      wait_for(g::q0);

      for (int i = 0; i < n; ++i) {
         buf[i] = group::grplist[i];
      }
      darray::copyin(g::q0, n, st.grplist, buf.data());
      wait_for(g::q0);

      for (int i = 0; i <= st.ngrp; ++i) {
         int j = 2 * i;
         buf[j] = group::igrp[j] - 1;
         buf[j + 1] = group::igrp[j + 1];
      }
      darray::copyin(g::q0, st.ngrp + 1, st.igrp, buf.data());
      darray::copyin(g::q0, st.ngrp + 1, st.grpmass, group::grpmass);
      wait_for(g::q0);

      std::vector<real> wgrpv((1 + st.ngrp) * (1 + st.ngrp));
      for (int i = 0; i <= st.ngrp; ++i) {
         for (int j = i; j <= st.ngrp; ++j) {
            real wg = group::wgrp[j + i * (1 + sizes::maxgrp)];
            wgrpv[j + i * (1 + st.ngrp)] = wg;
            wgrpv[i + j * (1 + st.ngrp)] = wg;
         }
      }
      darray::copyin(g::q0, wgrpv.size(), st.wgrp, wgrpv.data());
      wait_for(g::q0);
   }
}
}
