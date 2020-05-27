#pragma once
#include "energy_buffer.h"
#include "mdprec.h"


namespace tinker {
/**
 * \ingroup mdegv
 * \brief Kinetic energy on host.
 */
extern energy_prec eksum;
/**
 * \ingroup mdegv
 * \brief Kinetic energy tensor on host.
 */
extern energy_prec ekin[3][3];


//====================================================================//


/**
 * \ingroup mdegv
 * \brief
 * Zero out all of the counts, energies, gradients, and virials on device.
 */
void zero_egv(int vers);
/**
 * \ingroup mdegv
 * \brief
 * Zero out all of the counts, energies, gradients, and virials on device, with
 * parameter #rc_flag & calc::vmask.
 *
 * \see rc_flag
 * \see calc::vmask
 */
void zero_egv();


//====================================================================//


/**
 * \ingroup mdegv
 */
void scale_gradient(double scale, grad_prec* g0x, grad_prec* g0y,
                    grad_prec* g0z);
void sum_gradient(grad_prec* g0x, grad_prec* g0y, grad_prec* g0z,
                  const grad_prec* g1x, const grad_prec* g1y,
                  const grad_prec* g1z);
void sum_gradient(double scale, grad_prec* g0x, grad_prec* g0y, grad_prec* g0z,
                  const grad_prec* g1x, const grad_prec* g1y,
                  const grad_prec* g1z);
void scale_gradient_acc(double, grad_prec*, grad_prec*, grad_prec*);
void sum_gradient_acc(grad_prec*, grad_prec*, grad_prec*, const grad_prec*,
                      const grad_prec*, const grad_prec*);
void sum_gradient_acc(double, grad_prec*, grad_prec*, grad_prec*,
                      const grad_prec*, const grad_prec*, const grad_prec*);


//====================================================================//


/**
 * \ingroup mdegv
 * \brief
 * Copy total potential energy to another variable. Avoid accessing #esum
 * directly because it is not always available in the calculation, in which
 * case the output variable will not be changed.
 */
void copy_energy(int vers, energy_prec* eng);
/**
 * \ingroup mdegv
 * \brief
 * Copy the energy gradients from device to host.
 */
void copy_gradient(int vers, double* grdx, double* grdy, double* grdz,
                   const grad_prec* gx_src, const grad_prec* gy_src,
                   const grad_prec* gz_src);
/**
 * \ingroup mdegv
 * \brief
 * Copy the energy gradients from #gx, #gy, #gz to host.
 *
 * \see gx
 * \see gy
 * \see gz
 */
void copy_gradient(int vers, double* grdx, double* grdy, double* grdz);
/**
 * \ingroup mdegv
 * \brief
 * Copy virial tensor to another variable. Avoid accessing #vir directly
 * because it is not always available in the calculation, in which case the
 * output variable will not be changed.
 */
void copy_virial(int vers, virial_prec* virial);


//====================================================================//


void egv_data(rc_op op);
}
