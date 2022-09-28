#pragma once
#include "ff/precision.h"

namespace tinker {
/// \addtogroup pcg
/// \{

/// \page pcg Preconditioned Conjugate Gradient (PCG)
///
/// The PCG method in Tinker9 is implemented as follows.
///
/// 1. Goal: To solve Au = E, where A = 1/a + T.
///    - (a) E:   external field;
///    - (b) a:   polarizability;
///    - (c) -Tu: mutual induced field; i.e., \c ufield(u).
///
/// 2. Guess initial induced dipole: u.
///    - (a) Predict:      u = \c fromPrediction
///    - (b) Direct guess: u = aE
///    - (c) Zero:         u = 0
///
/// 3. Compute initial residual: r = E - Au
///    - (a) Predict:      r = (aE-u)/a + \c ufield(u)
///    - (b) Direct guess: r = \c ufield(aE)
///    - (c) Zero:         r = E
///
/// 4. p = Mr; s0 = rMr
///    - M is the preconditioner which approximates A.
///
/// 5. while not converged:
///    - (a) g = s0 / pAp, where Ap = a/p - \c ufield(p)
///    - (b) u = u + gp
///    - (c) r = r - gAp
///    - (d) s = rMr
///    - (e) b = s / s0
///    - (f) p = Mr + bp; s0 = s

/// \brief Calculates the direct induced dipoles \c 2b.
__global__
void pcgUdirV1(int n,                 ///< Number of atoms.
               const real* polarity,  ///< Polarizabilities.
               real (*udir)[3],       ///< [out] Direct induced dipoles \c 2b.
               const real (*field)[3] ///< External field.
);

/// \brief Calculates the direct induced dipoles \c 2b.
__global__
void pcgUdirV2(int n,                  ///< Number of atoms.
               const real* polarity,   ///< Polarizabilities.
               real (*udir)[3],        ///< [out] Direct induced d-dipoles \c 2b.
               real (*udirp)[3],       ///< [out] Direct induced p-dipoles \c 2b.
               const real (*field)[3], ///< External d-field.
               const real (*fieldp)[3] ///< External p-field.
);

/// \brief Calculates the initial residuals for the predictor \c 3a.
__global__
void pcgRsd0V1(int n,                    ///< Number of atoms.
               const real* polarity_inv, ///< Multiplicative inverse of polarizabilities.
               real (*rsd)[3],           ///< [out] Initial residuals.
               const real (*udir)[3],    ///< Direct induced dipoles \c 2b.
               const real (*uind)[3],    ///< Initial induced dipoles \c 2a.
               const real (*field)[3]    ///< Mutual field.
);

/// \brief Calculates the initial residuals for the predictor \c 3a.
__global__
void pcgRsd0V2(int n,                    ///< Number of atoms.
               const real* polarity_inv, ///< Multiplicative inverse of polarizabilities.
               real (*rsd)[3],           ///< [out] Initial d-residuals.
               real (*rsp)[3],           ///< [out] Initial p-residuals.
               const real (*udir)[3],    ///< Direct induced d-dipoles \c 2b.
               const real (*udip)[3],    ///< Direct induced p-dipoles \c 2b.
               const real (*uind)[3],    ///< Initial induced d-dipoles \c 2a.
               const real (*uinp)[3],    ///< Initial induced p-dipoles \c 2a.
               const real (*field)[3],   ///< Mutual d-field.
               const real (*fielp)[3]    ///< Mutual p-field.
);

/// \brief Calculates the initial residuals for the predictor \c 3a
/// using the <a href="https://doi.org/10.1021/acs.jpcb.2c04237">Exchange-Polarization Model.</a>
__global__
void pcgRsd0V3(int n,                       ///< Number of atoms.
               const real* polarity_inv,    ///< Multiplicative inverse of polarizabilities.
               real (*rsd)[3],              ///< [out] Initial residuals.
               const real (*udir)[3],       ///< Direct induced dipoles \c 2b.
               const real (*uind)[3],       ///< Initial induced dipoles \c 2a.
               const real (*field)[3],      ///< Mutual field.
               const real (*polscale)[3][3] ///< The \f$ I+kS^2 \f$ matrices for every atom (Eq. (7d)).
);

/// \}
}
