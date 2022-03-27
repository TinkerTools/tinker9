#include "ff/box.h"
#include "ff/energy.h"
#include "tinker9.h"
#include "tool/darray.h"
#include "tool/error.h"
#include <tinker/detail/atoms.hh>

namespace tinker {
void printError()
{
   Box p;
   boxGetDefault(p);
   boxSetTinkerModule(p);
   darray::copyout(g::q0, n, atoms::x, xpos);
   darray::copyout(g::q0, n, atoms::y, ypos);
   darray::copyout(g::q0, n, atoms::z, zpos);
   wait_for(g::q0);
   tinker_f_prterr();
}
}
