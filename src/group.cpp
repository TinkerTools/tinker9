#include "group.h"
#include "dev_array.h"
#include "md.h"
#include <tinker/detail/group.hh>
#include <tinker/detail/sizes.hh>


TINKER_NAMESPACE_BEGIN
void group_data(rc_op op)
{
   if (op & rc_dealloc) {
      auto& st = grp;
      device_array::deallocate(st.kgrp, st.grplist, st.igrp, st.grpmass,
                               st.wgrp);
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
      device_array::allocate(n, &st.kgrp, &st.grplist);
      device_array::allocate(1 + st.ngrp, &st.igrp, &st.grpmass);
      device_array::allocate((1 + st.ngrp) * (1 + st.ngrp), &st.wgrp);
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
      device_array::copyin(n, st.kgrp, buf.data());

      for (int i = 0; i < n; ++i) {
         buf[i] = group::grplist[i];
      }
      device_array::copyin(n, st.grplist, buf.data());

      for (int i = 0; i <= st.ngrp; ++i) {
         int j = 2 * i;
         buf[j] = group::igrp[j] - 1;
         buf[j + 1] = group::igrp[j + 1];
      }
      device_array::copyin(st.ngrp + 1, st.igrp, buf.data());

      device_array::copyin(st.ngrp + 1, st.grpmass, group::grpmass);

      std::vector<real> wgrpv((1 + st.ngrp) * (1 + st.ngrp));
      for (int i = 0; i <= st.ngrp; ++i) {
         for (int j = i; j <= st.ngrp; ++j) {
            real wg = group::wgrp[j + i * (1 + sizes::maxgrp)];
            wgrpv[j + i * (1 + st.ngrp)] = wg;
            wgrpv[i + j * (1 + st.ngrp)] = wg;
         }
      }
      device_array::copyin(wgrpv.size(), st.wgrp, wgrpv.data());
   }
}
TINKER_NAMESPACE_END
