#include "gpu/decl_dataop.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_pme.h"
#include "rc_cudart.h"
#include "util/math.h"
#include <ext/tinker/tinker_rt.h>
#include <vector>

TINKER_NAMESPACE_BEGIN
namespace gpu {
static void pme_op_destroy_(pme_st& st, pme_st*& dptr) {
  check_cudart(cudaFree(st.igrid));
  check_cudart(cudaFree(st.bsmod1));
  check_cudart(cudaFree(st.bsmod2));
  check_cudart(cudaFree(st.bsmod3));
  check_cudart(cudaFree(st.bsbuild));
  check_cudart(cudaFree(st.thetai1));
  check_cudart(cudaFree(st.thetai2));
  check_cudart(cudaFree(st.thetai3));
  check_cudart(cudaFree(st.qgrid));
  check_cudart(cudaFree(st.qfac));

  check_cudart(cudaFree(dptr));
}

static void pme_op_create_(pme_st& st, pme_st*& dptr, int nfft1, int nfft2,
                           int nfft3, int bsorder) {
  const size_t rs = sizeof(int);
  size_t size;

  st.nfft1 = nfft1;
  st.nfft2 = nfft2;
  st.nfft3 = nfft3;
  st.bsorder = bsorder;

  size = 3 * n * sizeof(int);
  check_cudart(cudaMalloc(&st.igrid, size));

  // see also subroutine moduli in pmestuf.f
  check_cudart(cudaMalloc(&st.bsmod1, rs * nfft1));
  check_cudart(cudaMalloc(&st.bsmod2, rs * nfft2));
  check_cudart(cudaMalloc(&st.bsmod3, rs * nfft3));
  int maxfft = max_of(nfft1, nfft2, nfft3);
  std::vector<double> array(bsorder);
  std::vector<double> bsarray(maxfft);
  double x = 0;
  TINKER_RT(bspline)(&x, &bsorder, array.data());
  for (int i = 0; i < maxfft; ++i) {
    bsarray[i] = 0;
  }
  assert(bsorder + 1 <= maxfft);
  for (int i = 0; i < bsorder; ++i) {
    bsarray[i + 1] = array[i];
  }
  std::vector<double> bsmodbuf(maxfft);
  TINKER_RT(dftmod)(bsmodbuf.data(), bsarray.data(), &nfft1, &bsorder);
  copyin_data(st.bsmod1, bsmodbuf.data(), nfft1);
  TINKER_RT(dftmod)(bsmodbuf.data(), bsarray.data(), &nfft2, &bsorder);
  copyin_data(st.bsmod2, bsmodbuf.data(), nfft2);
  TINKER_RT(dftmod)(bsmodbuf.data(), bsarray.data(), &nfft3, &bsorder);
  copyin_data(st.bsmod3, bsmodbuf.data(), nfft3);

  check_cudart(cudaMalloc(&st.bsbuild, rs * bsorder * bsorder));
  size = 4 * n * rs;
  check_cudart(cudaMalloc(&st.thetai1, size));
  check_cudart(cudaMalloc(&st.thetai2, size));
  check_cudart(cudaMalloc(&st.thetai3, size));
  size = nfft1 * nfft2 * nfft3 * rs;
  check_cudart(cudaMalloc(&st.qgrid, 2 * size));
  check_cudart(cudaMalloc(&st.qfac, size));

  size = sizeof(pme_st);
  check_cudart(cudaMalloc(&dptr, size));
  check_cudart(cudaMemcpy(dptr, &st, size, cudaMemcpyHostToDevice));
}

namespace detail_ {
std::vector<pme_st>& pme_objs() {
  static std::vector<pme_st> objs;
  return objs;
}

std::vector<pme_st*>& pme_deviceptrs() {
  static std::vector<pme_st*> ptrs;
  return ptrs;
}
}

int pme_open_unit(int nfft1, int nfft2, int nfft3, int bsorder) {
  bool found = false;
  int idx;
  for (idx = 0; found == false && idx < detail_::pme_objs().size(); ++idx) {
    auto& st = detail_::pme_objs()[idx];
    found = (nfft1 == st.nfft1 && nfft2 == st.nfft2 && nfft3 == st.nfft3 &&
             bsorder == st.bsorder);
  }

  if (!found) {
    detail_::pme_objs().emplace_back(pme_st());
    detail_::pme_deviceptrs().emplace_back(nullptr);
    auto& st = detail_::pme_objs()[idx];
    auto& dptr = detail_::pme_deviceptrs()[idx];
    pme_op_create_(st, dptr, nfft1, nfft2, nfft3, bsorder);
  }

  assert(detail_::pme_objs().size() == detail_::pme_deviceptrs().size());
  return idx;
}

pme_st& pme_obj(int index) { return detail_::pme_objs()[index]; }

pme_st* pme_deviceptr(int index) { return detail_::pme_deviceptrs()[index]; }

int epme_unit; // electrostatic
int ppme_unit; // polarization
int dpme_unit; // dispersion
}
TINKER_NAMESPACE_END

#include "gpu/e_mpole.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
void pme_data(int op) {
  if (op == op_destroy) {
    assert(detail_::pme_objs().size() == detail_::pme_deviceptrs().size());
    int idx = 0;
    while (idx < detail_::pme_objs().size()) {
      auto& st = detail_::pme_objs()[idx];
      auto& dptr = detail_::pme_deviceptrs()[idx];
      pme_op_destroy_(st, dptr);
      ++idx;
    }
    detail_::pme_objs().clear();
    detail_::pme_deviceptrs().clear();
    assert(detail_::pme_objs().size() == 0);
    assert(detail_::pme_deviceptrs().size() == 0);
  }

  if (op == op_create) {
    // electrostatics
    epme_unit = -1;
    if (use_empole()) {
      int typ;
      std::string typ_str;

      get_empole_type(typ, typ_str);
      if (typ == elec_ewald)
        epme_unit =
            pme_open_unit(pme::nefft1, pme::nefft2, pme::nefft3, pme::bseorder);
    }

    // polarization
    ppme_unit = -1;
    if (false) {
      ppme_unit =
          pme_open_unit(pme::nefft1, pme::nefft2, pme::nefft3, pme::bsporder);
    }

    // dispersion
    dpme_unit = -1;
    if (false) {
      dpme_unit =
          pme_open_unit(pme::ndfft1, pme::ndfft2, pme::ndfft3, pme::bsdorder);
    }
  }
}
}
TINKER_NAMESPACE_END
