#include "molecule.h"
#include "md.h"
#include "tool/darray.h"
#include <tinker/detail/molcul.hh>

namespace tinker {
void molecule_data(RcOp op)
{
   if (op & rc_dealloc) {
      auto& st = molecule;
      darray::deallocate(st.imol, st.kmol, st.molecule, st.molmass);
   }

   if (op & rc_alloc) {
      auto& st = molecule;
      darray::allocate(n, &st.imol, &st.kmol, &st.molecule, &st.molmass);
   }

   if (op & rc_init) {
      auto& st = molecule;

      std::vector<int> buf(2 * n);
      st.nmol = molcul::nmol;
      for (int i = 0; i < st.nmol; ++i) {
         int j = 2 * i;
         buf[j] = molcul::imol[j] - 1;
         buf[j + 1] = molcul::imol[j + 1];
      }
      darray::copyin(g::q0, st.nmol, st.imol, buf.data());
      wait_for(g::q0);
      for (int i = 0; i < n; ++i) {
         buf[i] = molcul::kmol[i] - 1;
      }
      darray::copyin(g::q0, n, st.kmol, buf.data());
      wait_for(g::q0);
      for (int i = 0; i < n; ++i) {
         buf[i] = molcul::molcule[i] - 1;
      }
      darray::copyin(g::q0, n, st.molecule, buf.data());
      wait_for(g::q0);
      st.totmass = molcul::totmass;
      darray::copyin(g::q0, st.nmol, st.molmass, molcul::molmass);
      wait_for(g::q0);
   }
}
}
