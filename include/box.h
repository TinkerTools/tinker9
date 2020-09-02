#pragma once
#include "glob.box.h"
#include "tool/rc_man.h"


namespace tinker {
void box_extent(double new_extent);


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
 * \ingroup box
 * \brief Similar to Tinker `lattice` subroutine.
 *
 * See doc/box/box.pdf for details.
 */
void box_lattice(Box& p, BoxShape sh, double a, double b, double c,
                 double alpha_deg, double beta_deg, double gamma_deg);


#define TINKER_IMAGE_LVEC_PARAMS  real3 lvec1, real3 lvec2, real3 lvec3
#define TINKER_IMAGE_LVEC_ARGS    lvec1, lvec2, lvec3
#define TINKER_IMAGE_RECIP_PARAMS real3 recipa, real3 recipb, real3 recipc
#define TINKER_IMAGE_RECIP_ARGS   recipa, recipb, recipc
#define TINKER_IMAGE_PARAMS                                                    \
   BoxShape box_shape, TINKER_IMAGE_LVEC_PARAMS, TINKER_IMAGE_RECIP_PARAMS
#define TINKER_IMAGE_ARGS                                                      \
   box_shape, TINKER_IMAGE_LVEC_ARGS, TINKER_IMAGE_RECIP_ARGS


void box_data(rc_op);
void box_data_acc(rc_op);
void box_copyin_acc();


/**
 * \ingroup box
 * \brief Get the volume of the PBC box.
 * \note This function may calculate the volume on-the-fly if there is not an
 * internal variable to save the volume, so the calculation must be cheap.
 * \note The volume is undefined for the non-PBC situation.
 */
real volbox();
}
