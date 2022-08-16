#include "ff/atom.h"
#include "ff/box.h"
#include "ff/energybuffer.h"
#include "ff/nblist.h"
#include "tool/darray.h"
#include "tool/error.h"
#include "tool/externfunc.h"
#include "tool/gpucard.h"
#include <tinker/detail/atomid.hh>
#include <tinker/detail/atoms.hh>
#include <tinker/detail/bound.hh>
#include <tinker/detail/usage.hh>

#include <cassert>

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
#else
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
TINKER_FVOID2(acc1, cu1, copyPosToXyz);
void copyPosToXyz()
{
   if CONSTEXPR (sizeof(pos_prec) == sizeof(real))
      return;

   TINKER_FCALL2(acc1, cu1, copyPosToXyz);
}

void copyPosToXyz(bool refreshNBList)
{
   copyPosToXyz();
   if (refreshNBList)
      nblistRefresh();
}

TINKER_FVOID2(acc1, cu1, bounds);
void bounds()
{
   if (not bound::use_bounds)
      return;

   TINKER_FCALL2(acc1, cu1, bounds);
   copyPosToXyz();
}

inline namespace v1 {
enum
{
   DCD_HEADER = 0,
   DCD_TDELTA = 10,
   DCD_USEBOX = 11,
   DCD_CTRL_LEN = 21,

   DCD_TITLE_NCHARS = 80,

   DCD_AX = 0,
   DCD_COS_G = 1,
   DCD_BX = 2,
   DCD_COS_B = 3,
   DCD_COS_A = 4,
   DCD_CX = 5,
   DCD_XTAL_LEN = 6,
};

enum class Archive
{
   NONE = 0,
   XYZ = 1,
   DCD = 2,
};
}

static int dcdControl[DCD_CTRL_LEN] = {0};
static std::vector<float> dcdx, dcdy, dcdz;
static Archive archive = Archive::NONE;

static void dcdReadIntoBuffer(void* buffer, int nbyte, std::ifstream& ipt)
{
   int size1, size2;
   ipt.read((char*)&size1, sizeof(int));
   if (nbyte > 0) assert(nbyte == size1);
   ipt.read((char*)buffer, size1);
   ipt.read((char*)&size2, sizeof(int));
}

void readFrameOpen(const std::string& filename, std::ifstream& ipt)
{
   // get file format type by inspection of first character
   char a1;
   ipt.open(filename);
   ipt >> a1;
   auto arc = Archive::NONE;
   if (a1 == ' ')
      arc = Archive::XYZ;
   else if ('0' <= a1 and a1 <= '9')
      arc = Archive::XYZ;
   else
      arc = Archive::DCD;

   if (arc == Archive::DCD) {
      ipt.close();
      ipt.open(filename, std::ios::in | std::ios::binary);

      // read header info along with title and number of atoms
      dcdReadIntoBuffer(dcdControl, sizeof(int) * DCD_CTRL_LEN, ipt);

      int dcdTitleRecordLen;
      ipt.read((char*)&dcdTitleRecordLen, sizeof(int));
      std::vector<char> titlebuf;
      titlebuf.resize(dcdTitleRecordLen + sizeof(int));
      ipt.read(titlebuf.data(), dcdTitleRecordLen + sizeof(int));

      int dcdNAtom;
      dcdReadIntoBuffer(&dcdNAtom, sizeof(int), ipt);
      assert(n == dcdNAtom);
      dcdx.resize(n);
      dcdy.resize(n);
      dcdz.resize(n);
   }

   archive = arc;
}

static void readFrameDCD(std::ifstream& ipt)
{
   if (dcdControl[DCD_USEBOX]) {
      double dcdXtal[DCD_XTAL_LEN];
      dcdReadIntoBuffer(dcdXtal, sizeof(double) * DCD_XTAL_LEN, ipt);
      double ax = dcdXtal[DCD_AX], bx = dcdXtal[DCD_BX], cx = dcdXtal[DCD_CX];
      double al = 90., be = 90., ga = 90.;
      if (dcdXtal[DCD_COS_A] != 0.0) al = std::acos(dcdXtal[DCD_COS_A]) * radian;
      if (dcdXtal[DCD_COS_B] != 0.0) be = std::acos(dcdXtal[DCD_COS_B]) * radian;
      if (dcdXtal[DCD_COS_G] != 0.0) ga = std::acos(dcdXtal[DCD_COS_G]) * radian;
      Box p;
      boxLattice(p, box_shape, ax, bx, cx, al, be, ga);
      boxSetCurrent(p);
   }

   dcdReadIntoBuffer(dcdx.data(), sizeof(float) * n, ipt);
   dcdReadIntoBuffer(dcdy.data(), sizeof(float) * n, ipt);
   dcdReadIntoBuffer(dcdz.data(), sizeof(float) * n, ipt);
   for (int i = 0; i < n; ++i) {
      atoms::x[i] = dcdx[i];
      atoms::y[i] = dcdy[i];
      atoms::z[i] = dcdz[i];
   }
}

static void readFrameXYZ(std::ifstream& ipt)
{
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
}

void readFrameCopyinToXyz(std::ifstream& ipt, bool& done)
{
   if (!ipt) done = true;
   if (done) return;

   if (archive == Archive::XYZ)
      readFrameXYZ(ipt);
   else if (archive == Archive::DCD)
      readFrameDCD(ipt);

   xyzData(RcOp::INIT);

   if (!ipt.good() or ipt.peek() == EOF) done = true;
}

void readFrameClose(std::ifstream& ipt)
{
   ipt.close();
   dcdx.clear();
   dcdy.clear();
   dcdz.clear();
   archive = Archive::NONE;
}
}
