#include "gpu/acc.h"
#include "gpu/decl_mdstate.h"
#include "gpu/e_polar.h"
#include "gpu/rc.h"
#include "util/format_print.h"
#include <ext/tinker/tinker_mod.h>
#include <ext/tinker/tinker_rt.h>

TINKER_NAMESPACE_BEGIN
namespace gpu {
// similar to uscale0a/uscale0b routines
static inline void diag_precond(const real (*rsd)[3], const real (*rsdp)[3],
                                real (*zrsd)[3], real (*zrsdp)[3]) {
  #pragma acc parallel loop independent\
              deviceptr(polarity,rsd,rsdp,zrsd,zrsdp)
  for (int i = 0; i < n; ++i) {
    real poli = polarity[i];
    #pragma acc loop independent
    for (int j = 0; j < 3; ++j) {
      zrsd[i][j] = poli * rsd[i][j];
      zrsdp[i][j] = poli * rsdp[i][j];
    }
  }
}

/**
 * PCG
 *
 * M = preconditioner
 * T = inv_alpha + Tu, -Tu = ufield
 *
 * r(0) = E - T u(0)
 * p(0) (or conj(0)) = M r(0)
 * ------------------------------
 * gamma (or a) = r(i) M r(i) / p(i) T p(i)
 * u(i+1) = u(i) + a p(i)
 * r(i+1) = r(i) - a T p(i)
 * beta (or b(i+1)) = r(i+1) M r(i+1) / r(i) M r(i)
 * p(i+1) = M r(r+1) + b(i+1) p(i)
 * ------------------------------
 *
 * subroutine induce0a in induce.f
 * rsd = r
 * zrsd = M r
 * conj = p
 * vec = T P
 */
void induce_mutual_pcg1(real* gpu_ud, real* gpu_up) {

  const int n3 = 3 * n;

  real(*uind)[3] = reinterpret_cast<real(*)[3]>(gpu_ud);
  real(*uinp)[3] = reinterpret_cast<real(*)[3]>(gpu_up);

  real(*field)[3] = work01_;
  real(*fieldp)[3] = work02_;
  real(*rsd)[3] = work03_;
  real(*rsdp)[3] = work04_;
  real(*zrsd)[3] = work05_;
  real(*zrsdp)[3] = work06_;
  real(*conj)[3] = work07_;
  real(*conjp)[3] = work08_;
  real(*vec)[3] = work09_;
  real(*vecp)[3] = work10_;

  const bool dirguess = polpcg::pcgguess;

  // zero out the induced dipoles at each site

  zero_array(&uind[0][0], n3);
  zero_array(&uinp[0][0], n3);

  // get the electrostatic field due to permanent multipoles

  dfield(&field[0][0], &fieldp[0][0]);

  // direct induced dipoles

  #pragma acc parallel loop independent\
              deviceptr(polarity,udir,udirp,field,fieldp)
  for (int i = 0; i < n; ++i) {
    real poli = polarity[i];
    #pragma acc loop independent
    for (int j = 0; j < 3; ++j) {
      udir[i][j] = poli * field[i][j];
      udirp[i][j] = poli * fieldp[i][j];
    }
  }

  if (dirguess) {
    copy_array(&uind[0][0], &udir[0][0], n3);
    copy_array(&uinp[0][0], &udirp[0][0], n3);
  }

  // initial residual r(0)

  // if do not use pcgguess, r(0) = E - T Zero = E
  // if use pcgguess, r(0) = E - (inv_alpha + Tu) alpha E
  //                       = E - E -Tu udir
  //                       = -Tu udir

  if (dirguess) {
    ufield(&udir[0][0], &udirp[0][0], &rsd[0][0], &rsdp[0][0]);
  } else {
    copy_array(&rsd[0][0], &field[0][0], n3);
    copy_array(&rsdp[0][0], &fieldp[0][0], n3);
  }

  // initial M r(0) and p(0)

  diag_precond(rsd, rsdp, zrsd, zrsdp);
  copy_array(&conj[0][0], &zrsd[0][0], n3);
  copy_array(&conjp[0][0], &zrsdp[0][0], n3);

  // initial r(0) M r(0)

  real sum, sump;
  sum = dotprod(&rsd[0][0], &zrsd[0][0], n3);
  sump = dotprod(&rsdp[0][0], &zrsdp[0][0], n3);

  // conjugate gradient iteration of the mutual induced dipoles

  const bool debug = inform::debug;
  const int politer = polpot::politer;
  const real poleps = polpot::poleps;
  const real debye = units::debye;
  const real pcgpeek = polpcg::pcgpeek;
  const int maxiter = 100; // see also subroutine induce0a in induce.f

  bool done = false;
  int iter = 0;
  real eps = 100;
  real epsold;

  while (!done) {
    ++iter;

    // vec = (inv_alpha + Tu) conj, field = -Tu conj
    // vec = inv_alpha * conj - field
    ufield(&conj[0][0], &conjp[0][0], &field[0][0], &fieldp[0][0]);
    #pragma acc parallel loop independent\
                deviceptr(polarity_inv,vec,vecp,conj,conjp,field,fieldp)
    for (int i = 0; i < n; ++i) {
      real poli_inv = polarity_inv[i];
      #pragma acc loop independent
      for (int j = 0; j < 3; ++j) {
        vec[i][j] = poli_inv * conj[i][j] - field[i][j];
        vecp[i][j] = poli_inv * conjp[i][j] - fieldp[i][j];
      }
    }

    // a <- p T p
    real a, ap;
    a = dotprod(&conj[0][0], &vec[0][0], n3);
    ap = dotprod(&conjp[0][0], &vecp[0][0], n3);
    // a <- r M r / p T p
    if (a != 0)
      a = sum / a;
    if (ap != 0)
      ap = sump / ap;

    // u <- u + a p
    // r <- r - a T p
    #pragma acc parallel loop independent\
                deviceptr(uind,uinp,conj,conjp,rsd,rsdp,vec,vecp)
    for (int i = 0; i < n; ++i) {
      #pragma acc loop independent
      for (int j = 0; j < 3; ++j) {
        uind[i][j] += a * conj[i][j];
        uinp[i][j] += ap * conjp[i][j];
        rsd[i][j] -= a * vec[i][j];
        rsdp[i][j] -= ap * vecp[i][j];
      }
    }

    // calculate/update M r
    diag_precond(rsd, rsdp, zrsd, zrsdp);

    real b, bp;
    real sum1, sump1;
    sum1 = dotprod(&rsd[0][0], &zrsd[0][0], n3);
    sump1 = dotprod(&rsdp[0][0], &zrsdp[0][0], n3);
    if (sum != 0)
      b = sum1 / sum;
    if (sump != 0)
      bp = sump1 / sump;
    sum = sum1;
    sump = sump1;

    // calculate/update p
    #pragma acc parallel loop independent\
             deviceptr(conj,conjp,zrsd,zrsdp)
    for (int i = 0; i < n; ++i) {
      #pragma acc loop independent
      for (int j = 0; j < 3; ++j) {
        conj[i][j] = zrsd[i][j] + b * conj[i][j];
        conjp[i][j] = zrsdp[i][j] + bp * conjp[i][j];
      }
    }

    real epsd;
    real epsp;
    epsd = dotprod(&rsd[0][0], &rsd[0][0], n3);
    epsp = dotprod(&rsdp[0][0], &rsdp[0][0], n3);

    epsold = eps;
    eps = REAL_MAX(epsd, epsp);
    eps = debye * REAL_SQRT(eps / n);

    if (debug) {
      if (iter == 1) {
        print(stdout,
              "\n Determination of SCF Induced Dipole Moments "
              ":{0:4s}Iter{0:4s}RMS Residual (Debye)\n",
              "");
      }
      print(stdout, "{0:>8d}{2:7s}{1:>16.10f}\n", iter, eps, "");
    }

    if (eps < poleps)
      done = true;
    if (eps > epsold)
      done = true;
    if (iter >= politer)
      done = true;

    // apply a "peek" iteration to the mutual induced dipoles

    if (done) {
      #pragma acc parallel loop independent\
                  deviceptr(polarity,uind,uinp,rsd,rsdp)
      for (int i = 0; i < n; ++i) {
        real term = pcgpeek * polarity[i];
        #pragma acc loop independent
        for (int j = 0; j < 3; ++j) {
          uind[i][j] += term * rsd[i][j];
          uinp[i][j] += term * rsdp[i][j];
        }
      }
    }
  }

  // print the results from the conjugate gradient iteration

  if (debug) {
    print(stdout,
          " Induced Dipoles :{2:4s}Iterations{0:>5d}{2:6s}RMS "
          "Residual{1:>15.10f}\n",
          iter, eps, "");
  }

  // terminate the calculation if dipoles failed to converge

  if (iter >= maxiter || eps > epsold) {
    TINKER_RT(prterr)();
    m_tinker_throw(" INDUCE  --  Warning, Induced Dipoles are not Converged");
  }
}
}
TINKER_NAMESPACE_END
