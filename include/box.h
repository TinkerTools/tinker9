#pragma once
#include "rc_man.h"
#include "realn.h"


namespace tinker {
/**
 * \ingroup box
 */
enum class BoxShape
{
   UNBOUND,
   ORTHO,
   MONO,
   TRI,
   OCT
};
constexpr BoxShape UNBOUND_BOX = BoxShape::UNBOUND;
constexpr BoxShape ORTHO_BOX = BoxShape::ORTHO;
constexpr BoxShape MONO_BOX = BoxShape::MONO;
constexpr BoxShape TRI_BOX = BoxShape::TRI;
constexpr BoxShape OCT_BOX = BoxShape::OCT;
TINKER_EXTERN BoxShape box_shape;
TINKER_EXTERN real3 lvec1, lvec2, lvec3;
TINKER_EXTERN real3 recipa, recipb, recipc;
void box_extent(double new_extent);


/**
 * \ingroup box
 * \brief Periodic boundary conditions (PBC).
 *
 * PBC in Tinker is defined by lengths of three axes and three angles:
 * `a-axis`, `b-axis`, `c-axis`, `alpha`, `beta`, and `gamma`.
 *
 * In cartesian coordinates:
 *    - `a-axis` is always aligned with the x axis, so the vector of `a-axis`
 *       is `(ax,ay,az) = (a,0,0)`;
 *    - `b-axis` does not have z component, so the vector of `b-axis` is
 *       `(bx,by,bz) = (bx,by,0)`;
 *    - `c-axis` is `(cx,cy,cz)`.
 *
 * Tinker has another 3x3 matrix `lvec`:
 *    - `lvec(:,1)` (or `lvec[0][]`): `(ax,bx,cx)`;
 *    - `lvec(:,2)` (or `lvec[1][]`): `(ay,by,cy)`;
 *    - `lvec(:,3)` (or `lvec[2][]`): `(az,bz,cz)`.
 *
 * Triclinic:
 *    - `lvec(:,1)` (or `lvec[0][]`): `(ax,bx,cx)`;
 *    - `lvec(:,2)` (or `lvec[1][]`): `( 0,by,cy)`;
 *    - `lvec(:,3)` (or `lvec[2][]`): `( 0, 0,cz)`.
 *
 * Monoclinic:
 *    - `lvec(:,1)` (or `lvec[0][]`): `(a,0,cx)`;
 *    - `lvec(:,2)` (or `lvec[1][]`): `(0,b, 0)`;
 *    - `lvec(:,3)` (or `lvec[2][]`): `(0,0,cz)`;
 *    - `alpha` and `gamma` are 90 degrees.
 *
 * Orthogonal:
 *    - `lvec(:,1)` (or `lvec[0][]`): `(a,0,0)`;
 *    - `lvec(:,2)` (or `lvec[1][]`): `(0,b,0)`;
 *    - `lvec(:,3)` (or `lvec[2][]`): `(0,0,c)`;
 *    - `alpha`, `beta`, and `gamma` are 90 degrees.
 *
 * Reciprocal lattice (`recip`), is the inverse of `lvec`:
 *    - Fortran: `recip(:,1)`, `recip(:,2)`, `recip(:,3)`;
 *    - C++: `recip[0][]`, `recip[1][]`, `recip[2][]`;
 *    - Fractional coordinates `fi (i=1,2,3) = recip(:,i) (xr,yr,zr)`;
 *    - Cartesian coordinates `wr (w=1,2,3 or x,y,z) = inv_recip(:,w) (f1,f2,f3)
         = lvec(:,w) (f1,f2,f3)`.
 */
struct Box
{
   BoxShape box_shape;
   real3 lvec1, lvec2, lvec3;
   real3 recipa, recipb, recipc;
};
void set_default_box(const Box&);
void get_default_box(Box&);
void set_recip_box(real3&, real3&, real3&, BoxShape, const real3&, const real3&,
                   const real3&);
void set_default_recip_box();
void get_box_axes_angles(const Box&, double& a, double& b, double& c,
                         double& alpha, double& beta, double& gamma);
void set_tinker_box_module(const Box&);
void get_tinker_box_module(Box&);
/**
 * \ingroup
 * \brief Similar to Tinker `lattice` subroutine.
 *
 * See doc/box/box.pdf for details.
 */
void box_lattice(Box& p, BoxShape sh, double a, double b, double c,
                 double alpha_deg, double beta_deg, double gamma_deg);


/**
 * \ingroup box
 * \brief Host pointer to the PBC boxes of a trajectory.
 */
TINKER_EXTERN Box* trajbox;


#define TINKER_IMAGE_LVEC_PARAMS  real3 lvec1, real3 lvec2, real3 lvec3
#define TINKER_IMAGE_LVEC_ARGS    lvec1, lvec2, lvec3
#define TINKER_IMAGE_RECIP_PARAMS real3 recipa, real3 recipb, real3 recipc
#define TINKER_IMAGE_RECIP_ARGS   recipa, recipb, recipc
#if 1
#   define TINKER_IMAGE_PARAMS                                                 \
      BoxShape box_shape, TINKER_IMAGE_LVEC_PARAMS, TINKER_IMAGE_RECIP_PARAMS
#   define TINKER_IMAGE_ARGS                                                   \
      box_shape, TINKER_IMAGE_LVEC_ARGS, TINKER_IMAGE_RECIP_ARGS
#else
#   define TINKER_IMAGE_PARAMS                                                 \
      TINKER_IMAGE_LVEC_PARAMS, TINKER_IMAGE_RECIP_PARAMS
#   define TINKER_IMAGE_ARGS TINKER_IMAGE_LVEC_ARGS, TINKER_IMAGE_RECIP_ARGS
#endif

void box_data(rc_op);


/**
 * \ingroup box
 * \brief Get the volume of the PBC box.
 * \note This function may calculate the volume on-the-fly if there is not an
 * internal variable to save the volume, so the calculation must be cheap.
 * \note The volume is undefined for the non-PBC situation.
 */
real volbox();
}
