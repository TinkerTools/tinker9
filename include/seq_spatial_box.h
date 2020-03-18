#pragma once
#include "macro.h"
#include "mathfunc.h"
#include "seq_def.h"


TINKER_NAMESPACE_BEGIN
SEQ_ROUTINE
inline void box_to_ixyz(int& restrict ix, int& restrict iy, int& restrict iz,
                        int px, int py, int pz, int boxid)
{
   ix = boxid >> (pz + py);         // ix = boxid / (dimz*dimy)
   boxid &= ((1 << (pz + py)) - 1); // boxid = boxid % (dimz*dimy)
   iy = boxid >> pz;                // iy = boxid / dimz
   iz = boxid & ((1 << pz) - 1);    // iz = boxid % dimz
}


SEQ_ROUTINE
inline int ixyz_to_box(int px, int py, int pz, int ix, int iy, int iz)
{
   int id = (ix << (pz + py)) + (iy << pz) + iz;
   return id;
}


#if 0
namespace spatial_v2 {
SEQ_ROUTINE
inline void box_to_ixyz(int& restrict ix, int& restrict iy, int& restrict iz,
                        int px, int py, int pz, int boxid)
{
   ix = 0;
   iy = 0;
   iz = 0;
   for (int i = 0; i < pz; ++i) {
      int zid = (boxid >> (2 * i + 0)) & (1 << i);
      int yid = (boxid >> (2 * i + 1)) & (1 << i);
      int xid = (boxid >> (2 * i + 2)) & (1 << i);
      iz += zid;
      iy += yid;
      ix += xid;
   }
   for (int i = pz; i < py; ++i) {
      int yid = (boxid >> (2 * i + 0)) & (1 << i);
      int xid = (boxid >> (2 * i + 1)) & (1 << i);
      iy += yid;
      ix += xid;
   }
   for (int i = py; i < px; ++i) {
      int xid = (boxid >> (2 * i + 0)) & (1 << i);
      ix += xid;
   }
}


SEQ_ROUTINE
inline int ixyz_to_box(int px, int py, int pz, int ix, int iy, int iz)
{
   int id = 0;
   for (int i = 0; i < pz; ++i) {
      int zid = (iz & (1 << i)) << (2 * i + 0);
      int yid = (iy & (1 << i)) << (2 * i + 1);
      int xid = (ix & (1 << i)) << (2 * i + 2);
      id += (zid + yid + xid);
   }
   for (int i = pz; i < py; ++i) {
      int yid = (iy & (1 << i)) << (2 * i + 0);
      int xid = (ix & (1 << i)) << (2 * i + 1);
      id += (yid + xid);
   }
   for (int i = py; i < px; ++i) {
      int xid = (ix & (1 << i)) << (2 * i + 0);
      id += xid;
   }
   return id;
}
}
#endif


SEQ_ROUTINE
inline void frac_to_ixyz(int& restrict ix, int& restrict iy, int& restrict iz,
                         int px, int py, int pz, real fx, real fy, real fz)
{
   ix = fx * (1 << px) + (1 << px) / 2; // ix = (fx+half) * 2^px;
   iy = fy * (1 << py) + (1 << py) / 2;
   iz = fz * (1 << pz) + (1 << pz) / 2;
}
TINKER_NAMESPACE_END
