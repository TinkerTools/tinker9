#include "gpu/mdstate.h"
#include "tinker.mod.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
int use_data = 0;
int n = 0;

real* esum;
real* vir;
real* mass;

real *x, *y, *z;
real *vx, *vy, *vz;
real *ax, *ay, *az;
real *gx, *gy, *gz;

box_st* box;
couple_st couple_obj_;
couple_st* couple;

nblist_st vlist_obj_;
nblist_st* vlist;

nblist_st elist_obj_;
nblist_st* elist;

void n_data() { n = atoms::n; }

void xyz_data(int op) {
  if ((use_xyz & use_data) == 0)
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(x));
    check_cudart(cudaFree(y));
    check_cudart(cudaFree(z));
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&x, size));
    check_cudart(cudaMalloc(&y, size));
    check_cudart(cudaMalloc(&z, size));
    copyin_data_1(x, atoms::x, n);
    copyin_data_1(y, atoms::y, n);
    copyin_data_1(z, atoms::z, n);
  }
}

void vel_data(int op) {
  if ((use_vel & use_data) == 0)
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(vx));
    check_cudart(cudaFree(vy));
    check_cudart(cudaFree(vz));
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&vx, size));
    check_cudart(cudaMalloc(&vy, size));
    check_cudart(cudaMalloc(&vz, size));
    copyin_data_n(0, 3, vx, moldyn::v, n);
    copyin_data_n(1, 3, vy, moldyn::v, n);
    copyin_data_n(2, 3, vz, moldyn::v, n);
  }
}

void accel_data(int op) {
  if ((use_accel & use_data) == 0)
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(ax));
    check_cudart(cudaFree(ay));
    check_cudart(cudaFree(az));
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&ax, size));
    check_cudart(cudaMalloc(&ay, size));
    check_cudart(cudaMalloc(&az, size));
    copyin_data_n(0, 3, ax, moldyn::a, n);
    copyin_data_n(1, 3, ay, moldyn::a, n);
    copyin_data_n(2, 3, az, moldyn::a, n);
  }
}

void mass_data(int op) {
  if ((use_mass & use_data) == 0)
    return;

  if (op == op_destroy) {
    check_cudart(cudaFree(mass));
  }

  if (op == op_create) {
    int n = atoms::n;
    size_t size = sizeof(real) * n;
    check_cudart(cudaMalloc(&mass, size));
    copyin_data_1(mass, atomid::mass, n);
  }
}

void energy_data(int op) {
  if (use_mass & (use_energy + use_grad + use_virial) == 0)
    return;

  if (op == op_destroy) {
    if (use_energy & use_data) {
      check_cudart(cudaFree(esum));
    }

    if (use_grad & use_data) {
      check_cudart(cudaFree(gx));
      check_cudart(cudaFree(gy));
      check_cudart(cudaFree(gz));
    }

    if (use_virial & use_data) {
      check_cudart(cudaFree(vir));
    }
  }

  if (op == op_create) {
    const size_t rs = sizeof(real);
    size_t size = 0;
    if (use_energy & use_data) {
      size = rs;
      check_cudart(cudaMalloc(&esum, size));
      copyin_data_1(esum, &energi::esum, 1);
    }

    if (use_grad & use_data) {
      int n = atoms::n;
      size = rs * n;
      check_cudart(cudaMalloc(&gx, size));
      check_cudart(cudaMalloc(&gy, size));
      check_cudart(cudaMalloc(&gz, size));
      copyin_data_n(0, 3, gx, deriv::desum, n);
      copyin_data_n(1, 3, gy, deriv::desum, n);
      copyin_data_n(2, 3, gz, deriv::desum, n);
    }

    if (use_virial & use_data) {
      size = rs * 9;
      check_cudart(cudaMalloc(&vir, size));
      copyin_data_1(vir, &virial::vir[0][0], 9);
    }
  }
}

void box_data(int op) {
  if (op == op_destroy) {
    check_cudart(cudaFree(box));
  }

  if (op == op_create) {
    int shape = box_null;
    if (boxes::orthogonal)
      shape = box_ortho;
    else if (boxes::monoclinic)
      shape = box_mono;
    else if (boxes::triclinic)
      shape = box_tri;
    else if (boxes::octahedron)
      shape = box_oct;

    size_t size = sizeof(box_st);
    check_cudart(cudaMalloc(&box, size));

    copyin_data_1(&box->xbox, &boxes::xbox, 1);
    copyin_data_1(&box->ybox, &boxes::ybox, 1);
    copyin_data_1(&box->zbox, &boxes::zbox, 1);
    copyin_data_1(&box->alpha, &boxes::alpha, 1);
    copyin_data_1(&box->beta, &boxes::beta, 1);
    copyin_data_1(&box->gamma, &boxes::gamma, 1);
    copyin_data_1(&box->lvec[0][0], &boxes::lvec[0][0], 9);
    copyin_data_1(&box->recip[0][0], &boxes::recip[0][0], 9);
    copyin_data_1(&box->shape, &shape, 1);
  }
}

void couple_data(int op) {
  if (op == op_destroy) {
    check_cudart(cudaFree(couple_obj_.n12));
    check_cudart(cudaFree(couple_obj_.n13));
    check_cudart(cudaFree(couple_obj_.n14));
    check_cudart(cudaFree(couple_obj_.n15));
    check_cudart(cudaFree(couple_obj_.i12));
    check_cudart(cudaFree(couple_obj_.i13));
    check_cudart(cudaFree(couple_obj_.i14));
    check_cudart(cudaFree(couple_obj_.i15));
    check_cudart(cudaFree(couple));
  }

  if (op == op_create) {
    const size_t rs = sizeof(int);
    size_t size;

    size = n * rs;
    check_cudart(cudaMalloc(&couple_obj_.n12, size));
    check_cudart(cudaMalloc(&couple_obj_.n13, size));
    check_cudart(cudaMalloc(&couple_obj_.n14, size));
    check_cudart(cudaMalloc(&couple_obj_.n15, size));
    size = couple_st::maxn12 * n * rs;
    check_cudart(cudaMalloc(&couple_obj_.i12, size));
    size = couple_st::maxn13 * n * rs;
    check_cudart(cudaMalloc(&couple_obj_.i13, size));
    size = couple_st::maxn14 * n * rs;
    check_cudart(cudaMalloc(&couple_obj_.i14, size));
    size = couple_st::maxn15 * n * rs;
    check_cudart(cudaMalloc(&couple_obj_.i15, size));

    // see also attach.f
    const int maxn13 = 3 * couple_st::maxn12;
    const int maxn14 = 9 * couple_st::maxn12;
    const int maxn15 = 27 * couple_st::maxn12;
    std::vector<int> nbuf, ibuf;
    nbuf.resize(n);

    size = couple_st::maxn12 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n12[i];
      nbuf[i] = nn;
      int base = i * couple_st::maxn12;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i12[i][j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_data_1(couple_obj_.n12, nbuf.data(), n);
    copyin_data_1(&couple_obj_.i12[0][0], ibuf.data(), size);

    size = couple_st::maxn13 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n13[i];
      nbuf[i] = nn;
      int base = i * couple_st::maxn13;
      int bask = i * maxn13;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i13[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_data_1(couple_obj_.n13, nbuf.data(), n);
    copyin_data_1(&couple_obj_.i13[0][0], ibuf.data(), size);

    size = couple_st::maxn14 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n14[i];
      nbuf[i] = nn;
      int base = i * couple_st::maxn14;
      int bask = i * maxn14;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i14[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_data_1(couple_obj_.n14, nbuf.data(), n);
    copyin_data_1(&couple_obj_.i14[0][0], ibuf.data(), size);

    size = couple_st::maxn15 * n;
    ibuf.resize(size);
    for (int i = 0; i < n; ++i) {
      int nn = couple::n15[i];
      nbuf[i] = nn;
      int base = i * couple_st::maxn15;
      int bask = i * maxn15;
      for (int j = 0; j < nn; ++j) {
        int k = couple::i15[bask + j];
        ibuf[base + j] = k - 1;
      }
    }
    copyin_data_1(couple_obj_.n15, nbuf.data(), n);
    copyin_data_1(&couple_obj_.i15[0][0], ibuf.data(), size);

    size = sizeof(couple_st);
    check_cudart(cudaMalloc(&couple, size));
    check_cudart(
        cudaMemcpy(couple, &couple_obj_, size, cudaMemcpyHostToDevice));
  }
}

void nblist_data(int op) {}
}
TINKER_NAMESPACE_END
