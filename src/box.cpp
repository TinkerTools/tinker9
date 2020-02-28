#include "box.h"
#include "dev_array.h"
#include "mathfunc.h"
#include "md.h"
#include "tinker_rt.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/boxes.hh>


TINKER_NAMESPACE_BEGIN
void set_default_box(const Box& p)
{
   lvec1 = p.lvec1;
   lvec2 = p.lvec2;
   lvec3 = p.lvec3;
   recipa = p.recipa;
   recipb = p.recipb;
   recipc = p.recipc;
}


void get_default_box(Box& p)
{
   p.lvec1 = lvec1;
   p.lvec2 = lvec2;
   p.lvec3 = lvec3;
   p.recipa = recipa;
   p.recipb = recipb;
   p.recipc = recipc;
}


void set_recip_box(real3 lvec1, real3 lvec2, real3 lvec3, real3& recipa,
                   real3& recipb, real3& recipc)
{
   recipc.x = 0;
   recipc.y = 0;
   recipc.z = 1.0 / lvec3.z;

   recipb.x = 0;
   recipb.y = 1.0 / lvec2.y;
   recipb.z = -lvec2.z / (lvec2.y * lvec3.z);

   recipa.x = 1.0 / lvec1.x;
   recipa.y = -lvec1.y / (lvec1.x * lvec2.y);
   recipa.z = lvec1.y * lvec2.z - lvec1.z * lvec2.y;
   recipa.z /= (lvec1.x * lvec2.y * lvec3.z);
}


#define DOT3(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
void set_tinker_box_module(const Box& p)
{
   if (bound::use_bounds) {
      double ax[3] = {p.lvec1.x, p.lvec2.x, p.lvec3.x};
      double bx[3] = {p.lvec1.y, p.lvec2.y, p.lvec3.y};
      double cx[3] = {p.lvec1.z, p.lvec2.z, p.lvec3.z};


      double xbox = std::sqrt(DOT3(ax, ax));
      double ybox = std::sqrt(DOT3(bx, bx));
      double zbox = std::sqrt(DOT3(cx, cx));
      double cos_a = DOT3(bx, cx) / (ybox * zbox);
      double cos_b = DOT3(cx, ax) / (zbox * xbox);
      double cos_c = DOT3(ax, bx) / (xbox * ybox);
      double a_deg = radian_dp * std::acos(cos_a);
      double b_deg = radian_dp * std::acos(cos_b);
      double c_deg = radian_dp * std::acos(cos_c);


      boxes::xbox = xbox;
      boxes::ybox = ybox;
      boxes::zbox = zbox;
      boxes::alpha = a_deg;
      boxes::beta = b_deg;
      boxes::gamma = c_deg;
      TINKER_RT(lattice)();
   }
}


void get_tinker_box_module(Box& p)
{
   const auto& r = boxes::recip;
   const auto& l = boxes::lvec;
   p.recipa = make_real3(r[0][0], r[0][1], r[0][2]);
   p.recipb = make_real3(r[1][0], r[1][1], r[1][2]);
   p.recipc = make_real3(r[2][0], r[2][1], r[2][2]);
   p.lvec1 = make_real3(l[0][0], l[0][1], l[0][2]);
   p.lvec2 = make_real3(l[1][0], l[1][1], l[1][2]);
   p.lvec3 = make_real3(l[2][0], l[2][1], l[2][2]);
}


void box_data(rc_op op)
{
   if (op & rc_dealloc) {
      if (calc::traj & rc_flag) {
         std::free(trajbox);
      } else {
         trajbox = nullptr;
      }
   }


   if (op & rc_alloc) {
      if (calc::traj & rc_flag) {
         trajbox = (Box*)std::malloc(sizeof(Box) * trajn);
      }
   }


   if (op & rc_init) {
      box_shape = UNBOUND_BOX;
      if (boxes::orthogonal)
         box_shape = ORTHO_BOX;
      else if (boxes::monoclinic)
         box_shape = MONO_BOX;
      else if (boxes::triclinic)
         box_shape = TRI_BOX;
      else if (boxes::octahedron)
         box_shape = OCT_BOX;

      Box p;
      get_tinker_box_module(p);
      set_default_box(p);
   }
}


real volbox()
{
   return lvec1.x * lvec2.y * lvec3.z;
}
TINKER_NAMESPACE_END
