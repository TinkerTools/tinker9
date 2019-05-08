#include "gpu/e.vdw.h"
#include "gpu/acc.h"
#include "gpu/image.h"
#include "gpu/mdstate.h"
#include "gpu/nblist.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
static void evdw_reduce_xyz_() {
  #pragma acc data deviceptr(x,y,z,ired,kred,xred,yred,zred)
  #pragma acc parallel loop async(queue_nb)
  for (int i = 0; i < n; ++i) {
    int iv = ired[i];
    real rdn = kred[i];
    xred[i] = rdn * (x[i] - x[iv]) + x[iv];
    yred[i] = rdn * (y[i] - y[iv]) + y[iv];
    zred[i] = rdn * (z[i] - z[iv]) + z[iv];
  }
}
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_vlist_build() {
  m_tinker_using_namespace;

  if (gpu::use_vdw_list()) {
    gpu::evdw_reduce_xyz_();
    gpu::nblist_construct(gpu::vlist_obj_, gpu::vlst);
  }
}

void tinker_gpu_vlist_update() {
  m_tinker_using_namespace;
  if (gpu::use_vdw_list()) {
    gpu::evdw_reduce_xyz_();
    gpu::nblist_update(gpu::vlist_obj_, gpu::vlst);
  }
}
}

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE, int VDWTYP>
void evdw_tmpl() {
  constexpr int do_e = USE & use_energy;
  constexpr int do_g = USE & use_grad;
  constexpr int do_v = USE & use_virial;
  constexpr int do_a = USE & use_analyz;

  const int maxnlst = vlist_obj_.maxnlst;
  const real off2 = REAL_SQ(vlist_obj_.cutoff);

  static real* vscale = (real*)malloc(n * sizeof(real));
  // In order to use firstprivate, must assign values here.
  for (int i = 0; i < n; ++i) {
    vscale[i] = 1;
  }

  #pragma acc data deviceptr(x,y,z,gx,gy,gz,vir,box,couple,vlst,\
                             vscales,\
                             ired,kred,xred,yred,zred,\
                             ev)\
   copyin(vscale[0:n])
  {
    #pragma acc serial async(queue_nb)
    { *ev = 0; }

    #pragma acc parallel loop firstprivate(vscale[0:n])
    for (int i = 0; i < n; ++i) {
      const int n12i = couple->n12[i];
      const int n13i = couple->n13[i];
      const int n14i = couple->n14[i];
      const int n15i = couple->n15[i];
      #pragma acc loop independent
      for (int j = 0; j < n12i; ++j)
        vscale[couple->i12[i][j]] = vscales[0]; // v2scale
      #pragma acc loop independent
      for (int j = 0; j < n13i; ++j)
        vscale[couple->i13[i][j]] = vscales[1]; // v3scale
      #pragma acc loop independent
      for (int j = 0; j < n14i; ++j)
        vscale[couple->i14[i][j]] = vscales[2]; // v4scale
      #pragma acc loop independent
      for (int j = 0; j < n15i; ++j)
        vscale[couple->i15[i][j]] = vscales[3]; // v5scale

      int iv = ired[i];
      real redi = kred[i];
      real rediv = 1 - redi;
      int it = jvdw[i];
      real xi, yi, zi;
      xi = xred[i];
      yi = yred[i];
      zi = zred[i];

      int nvlsti = vlst->nlst[i];
      int base = i * maxnlst;
      #pragma acc loop independent
      for (int kk = 0; kk < nvlsti; ++kk) {
        int k = vlst->lst[base + kk];
        int kv = ired[k];
        int kt = jvdw[k];
        real xr, yr, zr;
        xr = xi - xred[k];
        yr = yi - yred[k];
        zr = zi - zred[k];

        image(xr, yr, zr, box);
        real rik2 = xr * xr + yr * yr + zr * zr;
        bool incl = rik2 <= off2;
        real e;
        if_constexpr(do_e) {
          e = (incl ? 1 : 0);
          e *= vscale[k];
          #pragma acc atomic update
          *ev += e;
        }
      }

      #pragma acc loop independent
      for (int j = 0; j < n12i; ++j)
        vscale[couple->i12[i][j]] = 1;
      #pragma acc loop independent
      for (int j = 0; j < n13i; ++j)
        vscale[couple->i13[i][j]] = 1;
      #pragma acc loop independent
      for (int j = 0; j < n14i; ++j)
        vscale[couple->i14[i][j]] = 1;
      #pragma acc loop independent
      for (int j = 0; j < n15i; ++j)
        vscale[couple->i15[i][j]] = 1;
    }
  }
}
}
TINKER_NAMESPACE_END

extern "C" {
m_tinker_using_namespace;
TINKER_NONBONDED_GEN(tinker_gpu_evdw_lj, gpu::evdw_tmpl, gpu::evdw_lj);
TINKER_NONBONDED_GEN(tinker_gpu_evdw_buck, gpu::evdw_tmpl, gpu::evdw_buck);
TINKER_NONBONDED_GEN(tinker_gpu_evdw_mm3hb, gpu::evdw_tmpl, gpu::evdw_mm3hb);
TINKER_NONBONDED_GEN(tinker_gpu_evdw_hal, gpu::evdw_tmpl, gpu::evdw_hal);
TINKER_NONBONDED_GEN(tinker_gpu_evdw_gauss, gpu::evdw_tmpl, gpu::evdw_gauss);

#define TINKER_GPU_EVDW_GEN_(ver)                                              \
  void tinker_gpu_evdw##ver() {                                                \
    if (gpu::vdwtyp == gpu::evdw_lj)                                           \
      tinker_gpu_evdw_lj##ver();                                               \
    else if (gpu::vdwtyp == gpu::evdw_buck)                                    \
      tinker_gpu_evdw_buck##ver();                                             \
    else if (gpu::vdwtyp == gpu::evdw_mm3hb)                                   \
      tinker_gpu_evdw_mm3hb##ver();                                            \
    else if (gpu::vdwtyp == gpu::evdw_hal)                                     \
      tinker_gpu_evdw_hal##ver();                                              \
    else if (gpu::vdwtyp == gpu::evdw_gauss)                                   \
      tinker_gpu_evdw_gauss##ver();                                            \
  }
TINKER_GPU_EVDW_GEN_(0);
TINKER_GPU_EVDW_GEN_(1);
TINKER_GPU_EVDW_GEN_(3);
TINKER_GPU_EVDW_GEN_(4);
TINKER_GPU_EVDW_GEN_(5);
TINKER_GPU_EVDW_GEN_(6);
#undef TINKER_GPU_EVDW_GEN_
}