#pragma once
#include "energy_buffer.h"
#include "mdprec.h"


/**
 * \defgroup md_egv  Energy, Energy Gradient, and Virial Tensor
 * \ingroup md
 */


TINKER_NAMESPACE_BEGIN
/**
 * \ingroup md_egv
 * \brief Total potential energy on host.
 */
extern energy_prec esum;
/**
 * \ingroup md_egv
 * \brief Kinetic energy on host.
 */
extern energy_prec eksum;
/**
 * \ingroup md_egv
 * \brief Kinetic energy tensor on host.
 */
extern energy_prec ekin[3][3];
/**
 * \ingroup md_egv
 * \brief Total potential energy buffer on device.
 */
extern energy_buffer eng_buf;


/**
 * \ingroup md_egv
 * \brief Total potential energy gradients on device.
 */
extern grad_prec *gx, *gy, *gz;


/**
 * \ingroup md_egv
 * \brief Total virial tensor on host.
 */
extern virial_prec vir[9];
/**
 * \ingroup md_egv
 * \brief Total virial tensor buffer on device.
 */
extern virial_buffer vir_buf;


//====================================================================//


/**
 * \ingroup md_egv
 * \brief
 * Zero out all of the counts, energies, gradients, and virials on device.
 */
void zero_egv(int vers);
/**
 * \ingroup md_egv
 * \brief
 * Zero out all of the counts, energies, gradients, and virials on device, with
 * parameter #rc_flag & calc::vmask.
 *
 * \see rc_flag
 * \see calc::vmask
 */
void zero_egv();


/**
 * \ingroup md_egv
 * \brief Zero out the non-default gradients.
 */
void zero_gradient(DMFlag flag, size_t nelem, real* gx, real* gy, real* gz);
/**
 * \ingroup md_egv
 * \brief Zero out the non-default deterministic gradients.
 */
void zero_gradient(DMFlag flag, size_t nelem, fixed* gx, fixed* gy, fixed* gz);
void zero_gradient_acc(DMFlag, size_t, real*, real*, real*);
void zero_gradient_acc(DMFlag, size_t, fixed*, fixed*, fixed*);


//====================================================================//


/**
 * \ingroup md_egv
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


// sum up potential energies and virials on device and save on host
/**
 * \ingroup md_egv
 * \brief Sum the potential energies, virial tensors from device buffers and
 * save the results to host variables.
 *
 * \note The interactions counts are not summed by this function.
 *
 * \see esum
 * \see vir
 */
void sum_energy(int vers);


//====================================================================//


/**
 * \ingroup md_egv
 * \brief
 * Copy total potential energy to another variable. Avoid accessing #esum
 * directly because it is not always available in the calculation, in which
 * case the output variable will not be changed.
 */
void copy_energy(int vers, energy_prec* eng);
/**
 * \ingroup md_egv
 * \brief
 * Copy the energy gradients from device to host.
 */
void copy_gradient(int vers, double* grdx, double* grdy, double* grdz,
                   const grad_prec* gx_src, const grad_prec* gy_src,
                   const grad_prec* gz_src);
/**
 * \ingroup md_egv
 * \brief
 * Copy the energy gradients from #gx, #gy, #gz to host.
 *
 * \see gx
 * \see gy
 * \see gz
 */
void copy_gradient(int vers, double* grdx, double* grdy, double* grdz);
/**
 * \ingroup md_egv
 * \brief
 * Copy virial tensor to another variable. Avoid accessing #vir directly
 * because it is not always available in the calculation, in which case the
 * output variable will not be changed.
 */
void copy_virial(int vers, virial_prec* virial);


//====================================================================//


void egv_data(rc_op op);
TINKER_NAMESPACE_END
