#include "ff/atom.h"
#include "ff/box.h"
#include "ff/energybuffer.h"
#include "ff/nblist.h"
#include "tool/darray.h"
#include "tool/externfunc.h"
#include "tool/gpucard.h"
#include <cassert>
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/bound.hh>
#include <tinker/detail/usage.hh>

namespace tinker {
void nData(RcOp op)
{
   if (op & RcOp::DEALLOC) {
      trajn = -1;
      n = 0;
      padded_n = 0;
      nelem_buffer = 0;
   }

   if (op & RcOp::ALLOC) {
      n = atoms::n;
      padded_n = (n + WARP_SIZE - 1) / WARP_SIZE;
      padded_n *= WARP_SIZE;

      if (calc::traj & rc_flag) {
         // trajn must have been initialized by this point
         assert(trajn >= 0);
      }

#if TINKER_CUDART
      nelem_buffer = gpuMaxNParallel(idevice);
      nelem_buffer = pow2Ge(nelem_buffer);
#elif TINKER_HOST
      nelem_buffer = 1;
#endif

      if (usage::nuse != n) {
         TINKER_THROW("All atoms must be active.");
      }
   }
}

void massData(RcOp op)
{
   if (not(calc::mass & rc_flag))
      return;

   if (op & RcOp::DEALLOC) {
      darray::deallocate(mass, massinv);
   }

   if (op & RcOp::ALLOC) {
      darray::allocate(n, &mass, &massinv);
   }

   if (op & RcOp::INIT) {
      std::vector<double> mbuf(n);
      for (int i = 0; i < n; ++i)
         mbuf[i] = 1 / atomid::mass[i];
      darray::copyin(g::q0, n, massinv, mbuf.data());
      darray::copyin(g::q0, n, mass, atomid::mass);
      waitFor(g::q0);
   }
}

void xyzData(RcOp op)
{
   if (not(calc::xyz & rc_flag))
      return;

   if (op & RcOp::DEALLOC) {
      if (calc::traj & rc_flag) {
         darray::deallocate(trajx, trajy, trajz);
         x = nullptr;
         y = nullptr;
         z = nullptr;
      } else {
         trajx = nullptr;
         trajy = nullptr;
         trajz = nullptr;
         darray::deallocate(x, y, z);
         if (sizeof(pos_prec) == sizeof(real)) {
            xpos = nullptr;
            ypos = nullptr;
            zpos = nullptr;
         } else {
            darray::deallocate(xpos, ypos, zpos);
         }
      }
   }

   if (op & RcOp::ALLOC) {
      if (calc::traj & rc_flag) {
         darray::allocate(n * trajn, &trajx, &trajy, &trajz);
         x = trajx;
         y = trajy;
         z = trajz;
      } else {
         darray::allocate(n, &x, &y, &z);
         if (sizeof(pos_prec) == sizeof(real)) {
            xpos = (pos_prec*)x;
            ypos = (pos_prec*)y;
            zpos = (pos_prec*)z;
         } else {
            darray::allocate(n, &xpos, &ypos, &zpos);
         }
      }
   }

   if (op & RcOp::INIT) {
      if (calc::traj & rc_flag) {
         darray::copyin(g::q0, n, x, atoms::x);
         darray::copyin(g::q0, n, y, atoms::y);
         darray::copyin(g::q0, n, z, atoms::z);
      } else {
         darray::copyin(g::q0, n, xpos, atoms::x);
         darray::copyin(g::q0, n, ypos, atoms::y);
         darray::copyin(g::q0, n, zpos, atoms::z);
         copyPosToXyz();
      }
   }
}
}

namespace tinker {
TINKER_F2EXTN(cu, 0, acc, 1, void, copyPosToXyz);
void copyPosToXyz()
{
   TINKER_F2PICK(cu, acc, copyPosToXyz);
}

void copyPosToXyz(bool refreshNBList)
{
   copyPosToXyz_acc();
   if (refreshNBList)
      nblistRefresh();
}

TINKER_F2EXTN(cu, 0, acc, 1, void, boundsP1);
static void boundsP1()
{
   TINKER_F2PICK(cu, acc, boundsP1);
}

void bounds()
{
   if (not bound::use_bounds)
      return;

   boundsP1();
   copyPosToXyz();
}

void readFrameCopyinToXyz(std::istream& ipt, int& done)
{
   if (done)
      return;

   if (ipt) {
      std::string line;
      std::getline(ipt, line); // n and title
      std::getline(ipt, line); // either box size or first atom
      // 18.643000   18.643000   18.643000   90.000000   90.000000   90.000000
      //  1  O      8.733783    7.084710   -0.688468     1     2     3
      double l1, l2, l3, a1, a2, a3;
      int matched = std::sscanf(line.data(), "%lf%lf%lf%lf%lf%lf", &l1, &l2, &l3, &a1, &a2, &a3);
      int row = 0;
      int index;
      char name[32];
      double xr, yr, zr;
      if (matched == 6) {
         Box p;
         boxLattice(p, box_shape, l1, l2, l3, a1, a2, a3);
         boxSetCurrent(p);
      } else {
         std::sscanf(line.data(), "%d%s%lf%lf%lf", &index, name, &xr, &yr, &zr);
         index -= 1;
         atoms::x[index] = xr;
         atoms::y[index] = yr;
         atoms::z[index] = zr;
         row = 1;
      }

      for (int ir = row; ir < n; ++ir) {
         std::getline(ipt, line);
         std::sscanf(line.data(), "%d%s%lf%lf%lf", &index, name, &xr, &yr, &zr);
         index -= 1;
         atoms::x[index] = xr;
         atoms::y[index] = yr;
         atoms::z[index] = zr;
      }

      xyzData(RcOp::INIT);
   }

   if (ipt.peek() == EOF)
      done = true;
}
}
