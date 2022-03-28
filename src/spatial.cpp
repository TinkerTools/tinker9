#include "ff/spatial.h"
#include "ff/box.h"
#include "ff/nblist.h"
#include "tool/darray.h"
#include "tool/error.h"

namespace tinker {
/// Order of "roll-cut": x-y-z-x-...
MAYBE_UNUSED
static void spatial_cut_v1(int& px, int& py, int& pz, int level)
{
   px = (level + 2) / 3;
   py = (level + 1) / 3;
   pz = (level + 0) / 3;
}
/// Order of "roll-cut": z-y-x-z-...
MAYBE_UNUSED
static void spatial_cut_v2(int& px, int& py, int& pz, int level)
{
   px = (level + 0) / 3;
   py = (level + 1) / 3;
   pz = (level + 2) / 3;
}
/// Order of "roll-cut": always cuts "the longest axis".
static void spatial_cut_v3(int& px, int& py, int& pz, int level)
{
   // triclinic frac(1,1,1) -> cart(x,y,z)
   // x = (fz * l1.z + fy * l1.y + fx * l1.x)
   // y = (fz * l2.z + fy * l2.y)
   // z = (fz * l3.z)

   double3 l1 = make_double3(lvec1.x, lvec1.y, lvec1.z);
   double3 l2 = make_double3(lvec2.x, lvec2.y, lvec2.z);
   double3 l3 = make_double3(lvec3.x, lvec3.y, lvec3.z);
   px = 0;
   py = 0;
   pz = 0;
   const double ratio = 0.95;
   for (int i = 0; i < level; ++i) {
      double xx = l1.z + l1.y + l1.x;
      double yy = l2.z + l2.y;
      double zz = l3.z;

      if ((zz > ratio * xx) && (zz > ratio * yy)) {
         // if z is approximately the longest, cut c-axis by half
         l1.z /= 2;
         l2.z /= 2;
         l3.z /= 2;
         pz += 1;
      } else if (yy > ratio * xx) {
         // if y is approximately the longest, cut b-axis by half
         l1.y /= 2;
         l2.y /= 2;
         l3.y /= 2;
         py += 1;
      } else {
         // if x is longest, cut a-axis by half
         l1.x /= 2;
         l2.x /= 2;
         l3.x /= 2;
         px += 1;
      }
   }
}
void spatial1_cut(int& px, int& py, int& pz, int level)
{
   spatial_cut_v3(px, py, pz, level);
}
}
