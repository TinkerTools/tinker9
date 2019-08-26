#include "acc_seq.h"
#include "array.h"
#include "e_polar.h"
#include "error.h"
#include "ext/tinker/detail/inform.hh"
#include "ext/tinker/detail/polpcg.hh"
#include "ext/tinker/detail/polpot.hh"
#include "ext/tinker/detail/units.hh"
#include "io_print.h"
#include "md.h"
#include "nblist.h"
#include "tinker_rt.h"

TINKER_NAMESPACE_BEGIN
// similar to uscale0a/uscale0b routines
// the preconditioner is the diagnoal matrix
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

// similar to uscale0a/uscale0b routines
// the preconditioner is the sparse diagnoal matrix
static inline void sparse_diag_precond_apply(const real (*rsd)[3],
                                             const real (*rsdp)[3],
                                             real (*zrsd)[3],
                                             real (*zrsdp)[3]) {
  #pragma acc parallel loop independent\
              deviceptr(polarity,rsd,rsdp,zrsd,zrsdp)
  for (int i = 0; i < n; ++i) {
    real poli = udiag * polarity[i];
    #pragma acc loop independent
    for (int j = 0; j < 3; ++j) {
      zrsd[i][j] = poli * rsd[i][j];
      zrsdp[i][j] = poli * rsdp[i][j];
    }
  }

  const int maxnlst = ulist_unit->maxnlst;
  const auto* ulst = ulist_unit.deviceptr();

  #pragma acc parallel loop independent\
              deviceptr(rsd,rsdp,zrsd,zrsdp,mindex,minv,ulst)
  for (int i = 0; i < n; ++i) {
    int m00 = mindex[i];
    int nulsti = ulst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop independent
    for (int kk = 0; kk < nulsti; ++kk) {
      int m = m00 + kk * 6;
      int k = ulst->lst[base + kk];
      real m0 = minv[m];
      real m1 = minv[m + 1];
      real m2 = minv[m + 2];
      real m3 = minv[m + 3];
      real m4 = minv[m + 4];
      real m5 = minv[m + 5];

      #pragma acc atomic update
      zrsd[i][0] += m0 * rsd[k][0] + m1 * rsd[k][1] + m2 * rsd[k][2];
      #pragma acc atomic update
      zrsd[i][1] += m1 * rsd[k][0] + m3 * rsd[k][1] + m4 * rsd[k][2];
      #pragma acc atomic update
      zrsd[i][2] += m2 * rsd[k][0] + m4 * rsd[k][1] + m5 * rsd[k][2];

      #pragma acc atomic update
      zrsd[k][0] += m0 * rsd[i][0] + m1 * rsd[i][1] + m2 * rsd[i][2];
      #pragma acc atomic update
      zrsd[k][1] += m1 * rsd[i][0] + m3 * rsd[i][1] + m4 * rsd[i][2];
      #pragma acc atomic update
      zrsd[k][2] += m2 * rsd[i][0] + m4 * rsd[i][1] + m5 * rsd[i][2];

      #pragma acc atomic update
      zrsdp[i][0] += m0 * rsdp[k][0] + m1 * rsdp[k][1] + m2 * rsdp[k][2];
      #pragma acc atomic update
      zrsdp[i][1] += m1 * rsdp[k][0] + m3 * rsdp[k][1] + m4 * rsdp[k][2];
      #pragma acc atomic update
      zrsdp[i][2] += m2 * rsdp[k][0] + m4 * rsdp[k][1] + m5 * rsdp[k][2];

      #pragma acc atomic update
      zrsdp[k][0] += m0 * rsdp[i][0] + m1 * rsdp[i][1] + m2 * rsdp[i][2];
      #pragma acc atomic update
      zrsdp[k][1] += m1 * rsdp[i][0] + m3 * rsdp[i][1] + m4 * rsdp[i][2];
      #pragma acc atomic update
      zrsdp[k][2] += m2 * rsdp[i][0] + m4 * rsdp[i][1] + m5 * rsdp[i][2];
    }
  }
}

static inline void sparse_diag_precond_build(const real (*rsd)[3],
                                             const real (*rsdp)[3],
                                             real (*zrsd)[3],
                                             real (*zrsdp)[3]) {
  const auto* nulst = ulist_unit->nlst;
  #pragma acc serial deviceptr(mindex, nulst)
  {
    int m = 0;
    for (int i = 0; i < n; ++i) {
      mindex[i] = m;
      m += 6 * nulst[i];
    }
  }

  const int maxnlst = ulist_unit->maxnlst;
  const auto* ulst = ulist_unit.deviceptr();

  const auto* polargroup = polargroup_unit.deviceptr();

  auto bufsize = EnergyBuffer::estimate_size(n);

  static std::vector<real> uscalebuf;
  uscalebuf.resize(n, 1);
  real* uscale = uscalebuf.data();

  #pragma acc parallel loop gang num_gangs(bufsize) independent\
              deviceptr(mindex,minv,ulst,box,\
              polargroup,x,y,z,polarity,pdamp,thole)\
              firstprivate(uscale[0:n])
  for (int i = 0; i < n; ++i) {

    // set exclusion coefficients for connected atoms

    const int np11i = polargroup->np11[i];
    const int np12i = polargroup->np12[i];
    const int np13i = polargroup->np13[i];
    const int np14i = polargroup->np14[i];

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      uscale[polargroup->ip11[i][j]] = u1scale;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      uscale[polargroup->ip12[i][j]] = u2scale;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      uscale[polargroup->ip13[i][j]] = u3scale;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      uscale[polargroup->ip14[i][j]] = u4scale;

    real xi = x[i];
    real yi = y[i];
    real zi = z[i];
    real poli = polarity[i];
    //
    real pdi = pdamp[i];
    real pti = thole[i];
    //

    int nulsti = ulst->nlst[i];
    int base = i * maxnlst;
    #pragma acc loop independent
    for (int kk = 0; kk < nulsti; ++kk) {
      int k = ulst->lst[base + kk];
      real xr = x[k] - xi;
      real yr = y[k] - yi;
      real zr = z[k] - zi;
      image(xr, yr, zr, box);
      real r2 = xr * xr + yr * yr + zr * zr;
      real r = REAL_SQRT(r2);

      //
      real scale3 = uscale[k];
      real scale5 = uscale[k];
      real damp = pdi * pdamp[k];
      if (damp != 0) {
        real pgamma = REAL_MIN(pti, thole[k]);
        damp = -pgamma * REAL_CUBE(r / damp);
        if (damp > -50) {
          real expdamp = REAL_EXP(damp);
          scale3 *= (1 - expdamp);
          scale5 *= (1 - expdamp * (1 - damp));
        }
      }
      //
      real polik = poli * polarity[k];
      real rr3 = scale3 * polik * REAL_RECIP(r * r2);
      real rr5 = 3 * scale5 * polik * REAL_RECIP(r * r2 * r2);

      int m = mindex[i] + 6 * kk;
      minv[m] = rr5 * xr * xr - rr3;
      minv[m + 1] = rr5 * xr * yr;
      minv[m + 2] = rr5 * xr * zr;
      minv[m + 3] = rr5 * yr * yr - rr3;
      minv[m + 4] = rr5 * yr * zr;
      minv[m + 5] = rr5 * zr * zr - rr3;
    } // end for (int kk)

    #pragma acc loop independent
    for (int j = 0; j < np11i; ++j)
      uscale[polargroup->ip11[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np12i; ++j)
      uscale[polargroup->ip12[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np13i; ++j)
      uscale[polargroup->ip13[i][j]] = 1;
    #pragma acc loop independent
    for (int j = 0; j < np14i; ++j)
      uscale[polargroup->ip14[i][j]] = 1;
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
  // use sparse matrix preconditioner
  // or just use diagonal matrix preconditioner
  const bool sparse_prec = polpcg::pcgprec;

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

  if (sparse_prec) {
    sparse_diag_precond_build(rsd, rsdp, zrsd, zrsdp);
    sparse_diag_precond_apply(rsd, rsdp, zrsd, zrsdp);
  } else {
    diag_precond(rsd, rsdp, zrsd, zrsdp);
  }
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
    if (sparse_prec)
      sparse_diag_precond_apply(rsd, rsdp, zrsd, zrsdp);
    else
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
              "\n Determination of SCF Induced Dipole Moments\n\n"
              "{0:4s}Iter{0:4s}RMS Residual (Debye)\n\n",
              "");
      }
      print(stdout, "{0:>8d}{2:8s}{1:<16.10f}\n", iter, eps, "");
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
    TINKER_THROW(" INDUCE  --  Warning, Induced Dipoles are not Converged");
  }
}
TINKER_NAMESPACE_END
