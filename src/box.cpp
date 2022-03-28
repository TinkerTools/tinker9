#include "ff/box.h"
#include "math/inc.h"
#include "md/inc.h"
#include "tinker9.h"
#include "tool/darray.h"
#include <tinker/detail/bound.hh>
#include <tinker/detail/boxes.hh>

namespace tinker {
void boxExtent(double new_extent)
{
   if (box_shape != BoxShape::UNBOUND)
      return;
   if (lvec1.x >= new_extent)
      return;

   double w1x = 1.0 / new_extent;
   lvec1.x = new_extent;
   lvec2.y = new_extent;
   lvec3.z = new_extent;
   recipa.x = w1x;
   recipb.y = w1x;
   recipc.z = w1x;
}

extern void boxCopyin_acc();
void boxCopyin()
{
   boxCopyin_acc();
}

void boxSetDefault(const Box& p)
{
   box_shape = p.box_shape;
   lvec1 = p.lvec1;
   lvec2 = p.lvec2;
   lvec3 = p.lvec3;
   recipa = p.recipa;
   recipb = p.recipb;
   recipc = p.recipc;

   boxCopyin();
}

void boxGetDefault(Box& p)
{
   p.box_shape = box_shape;
   p.lvec1 = lvec1;
   p.lvec2 = lvec2;
   p.lvec3 = lvec3;
   p.recipa = recipa;
   p.recipb = recipb;
   p.recipc = recipc;
}

static void boxSetRecip(real3& recipa, real3& recipb, real3& recipc, BoxShape box_shape,
   const real3& lvec1, const real3& lvec2, const real3& lvec3)
{
   if (box_shape == BoxShape::ORTHO) {
      recipc.x = 0;
      recipc.y = 0;
      recipc.z = 1.0 / lvec3.z;

      recipb.x = 0;
      recipb.y = 1.0 / lvec2.y;
      recipb.z = 0;

      recipa.x = 1.0 / lvec1.x;
      recipa.y = 0;
      recipa.z = 0;
   } else if (box_shape == BoxShape::MONO) {
      recipc.x = 0;
      recipc.y = 0;
      recipc.z = 1.0 / lvec3.z;

      recipb.x = 0;
      recipb.y = 1.0 / lvec2.y;
      recipb.z = 0;

      recipa.x = 1.0 / lvec1.x;
      recipa.y = 0;
      recipa.z = -lvec1.z / (lvec1.x * lvec3.z);
   } else if (box_shape == BoxShape::TRI) {
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
   } else if (box_shape == BoxShape::OCT) {
      recipc.x = 0;
      recipc.y = 0;
      recipc.z = 1.0 / lvec1.x;

      recipb.x = 0;
      recipb.y = 1.0 / lvec1.x;
      recipb.z = 0;

      recipa.x = 1.0 / lvec1.x;
      recipa.y = 0;
      recipa.z = 0;
   }
}

void boxSetDefaultRecip()
{
   boxSetRecip(recipa, recipb, recipc, box_shape, lvec1, lvec2, lvec3);

   boxCopyin();
}

void boxGetAxesAngles(
   const Box& p, double& a, double& b, double& c, double& alpha, double& beta, double& gamma)
{
   auto DOT3 = [](const double* a, const double* b) -> double {
      return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
   };

   double ax[3] = {p.lvec1.x, p.lvec2.x, p.lvec3.x};
   double bx[3] = {p.lvec1.y, p.lvec2.y, p.lvec3.y};
   double cx[3] = {p.lvec1.z, p.lvec2.z, p.lvec3.z};

   double xbox = std::sqrt(DOT3(ax, ax));
   double ybox = std::sqrt(DOT3(bx, bx));
   double zbox = std::sqrt(DOT3(cx, cx));
   double cos_a = DOT3(bx, cx) / (ybox * zbox);
   double cos_b = DOT3(cx, ax) / (zbox * xbox);
   double cos_c = DOT3(ax, bx) / (xbox * ybox);
   double a_deg = (180 / M_PI) * std::acos(cos_a);
   double b_deg = (180 / M_PI) * std::acos(cos_b);
   double c_deg = (180 / M_PI) * std::acos(cos_c);

   a = xbox;
   b = ybox;
   c = zbox;
   alpha = a_deg;
   beta = b_deg;
   gamma = c_deg;
}

void boxSetTinkerModule(const Box& p)
{
   if (p.box_shape == BoxShape::UNBOUND)
      return;

   boxes::orthogonal = 0;
   boxes::monoclinic = 0;
   boxes::triclinic = 0;
   boxes::octahedron = 0;
   if (box_shape == BoxShape::ORTHO)
      boxes::orthogonal = 1;
   else if (box_shape == BoxShape::MONO)
      boxes::monoclinic = 1;
   else if (box_shape == BoxShape::TRI)
      boxes::triclinic = 1;
   else if (box_shape == BoxShape::OCT)
      boxes::octahedron = 1;

   double xbox, ybox, zbox, a_deg, b_deg, c_deg;
   boxGetAxesAngles(p, xbox, ybox, zbox, a_deg, b_deg, c_deg);

   boxes::xbox = xbox;
   boxes::ybox = ybox;
   boxes::zbox = zbox;
   boxes::alpha = a_deg;
   boxes::beta = b_deg;
   boxes::gamma = c_deg;
   tinker_f_lattice();
}

static void boxGetTinkerModule(Box& p)
{
   if (!bound::use_bounds) {
      p.box_shape = BoxShape::UNBOUND;
      p.lvec1 = make_real3(0, 0, 0);
      p.lvec2 = make_real3(0, 0, 0);
      p.lvec3 = make_real3(0, 0, 0);
      p.recipa = make_real3(0, 0, 0);
      p.recipb = make_real3(0, 0, 0);
      p.recipc = make_real3(0, 0, 0);
      return;
   }

   if (boxes::orthogonal)
      p.box_shape = BoxShape::ORTHO;
   else if (boxes::monoclinic)
      p.box_shape = BoxShape::MONO;
   else if (boxes::triclinic)
      p.box_shape = BoxShape::TRI;
   else if (boxes::octahedron)
      p.box_shape = BoxShape::OCT;

   const auto& r = boxes::recip;
   const auto& l = boxes::lvec;
   p.recipa.x = r[0][0];
   p.recipa.y = r[0][1];
   p.recipa.z = r[0][2];
   p.recipb.x = 0; // r[1][0];
   p.recipb.y = r[1][1];
   p.recipb.z = r[1][2];
   p.recipc.x = 0; // r[2][0];
   p.recipc.y = 0; // r[2][1];
   p.recipc.z = r[2][2];
   p.lvec1.x = l[0][0];
   p.lvec1.y = l[0][1];
   p.lvec1.z = l[0][2];
   p.lvec2.x = 0; // l[1][0];
   p.lvec2.y = l[1][1];
   p.lvec2.z = l[1][2];
   p.lvec3.x = 0; // l[1][1];
   p.lvec3.y = 0; // l[1][2];
   p.lvec3.z = l[2][2];
}

void boxLattice(Box& p, BoxShape sh, double a, double b, double c, double alpha_deg,
   double beta_deg, double gamma_deg)
{
   p.box_shape = sh;
   if (sh == BoxShape::TRI) {
      double alpha = alpha_deg * M_PI / 180;
      double beta = beta_deg * M_PI / 180;
      double gamma = gamma_deg * M_PI / 180;
      double cos_alpha = std::cos(alpha);
      double cos_beta = std::cos(beta);
      double cos_gamma = std::cos(gamma);
      double sin_gamma = std::sin(gamma);
      double bx = b * cos_gamma;
      double by = b * sin_gamma;
      double cx = c * cos_beta;
      double cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
      double cz = c / sin_gamma *
         std::sqrt(2 * cos_alpha * cos_beta * cos_gamma + sin_gamma * sin_gamma -
            cos_alpha * cos_alpha - cos_beta * cos_beta);
      p.lvec1 = make_real3(a, bx, cx);
      p.lvec2 = make_real3(0, by, cy);
      p.lvec3 = make_real3(0, 0, cz);
   } else if (sh == BoxShape::MONO) {
      double beta = beta_deg * M_PI / 180;
      double cos_beta = std::cos(beta);
      double sin_beta = std::sin(beta);
      double cx = c * cos_beta;
      double cz = c * sin_beta;
      p.lvec1 = make_real3(a, 0, cx);
      p.lvec2 = make_real3(0, b, 0);
      p.lvec3 = make_real3(0, 0, cz);
   } else if (sh == BoxShape::ORTHO) {
      p.lvec1 = make_real3(a, 0, 0);
      p.lvec2 = make_real3(0, b, 0);
      p.lvec3 = make_real3(0, 0, c);
   } else if (sh == BoxShape::OCT) {
      p.lvec1 = make_real3(a, 0, 0);
      p.lvec2 = make_real3(0, a, 0);
      p.lvec3 = make_real3(0, 0, a);
   }
   boxSetRecip(p.recipa, p.recipb, p.recipc, p.box_shape, p.lvec1, p.lvec2, p.lvec3);
}

extern void boxData_acc(RcOp);
void boxData(RcOp op)
{
   boxData_acc(op);

   if (op & rc_dealloc) {
      if (calc::traj & rc_flag) {
         std::free(trajbox);
      } else {
         trajbox = nullptr;
      }
      box_shape = BoxShape::UNBOUND;
   }

   if (op & rc_alloc) {
      if (calc::traj & rc_flag) {
         trajbox = (Box*)std::malloc(sizeof(Box) * trajn);
      }
   }

   if (op & rc_init) {
      Box p;
      boxGetTinkerModule(p);
      boxSetDefault(p);
   }
}

real boxVolume()
{
   real ans = lvec1.x * lvec2.y * lvec3.z;
   if (box_shape == BoxShape::OCT)
      ans *= 0.5f;
   return ans;
}
}
