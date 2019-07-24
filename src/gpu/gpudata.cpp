#include "gpu/gpudata.h"
#include "gpu/gpu.h"

extern "C" {
void tinker_gpu_data_create() {
  m_tinker_using_namespace;
  using namespace gpu;
  const rc_t rc = static_cast<rc_t>(rc_alloc | rc_copyin);
  gpu_data(rc);
}

void tinker_gpu_data_create_() { tinker_gpu_data_create_(); }

void tinker_gpu_data_destroy() {
  m_tinker_using_namespace;
  using namespace gpu;
  const rc_t rc = rc_dealloc;
  gpu_data(rc);
}

void tinker_gpu_data_destroy_() { tinker_gpu_data_destroy(); }
}

TINKER_NAMESPACE_BEGIN
namespace gpu {
extern void n_data(rc_t);

extern void xyz_data(rc_t);
extern void vel_data(rc_t);
extern void mass_data(rc_t);

extern void potential_data(rc_t);

extern void box_data(rc_t);
extern void couple_data(rc_t);
extern void nblist_data(rc_t);

extern void random_data(rc_t);

extern void md_data(rc_t);

void gpu_data(rc_t rc) {
  rc_man<n_data> n42_{rc};

  rc_man<xyz_data> xyz42_{rc};
  rc_man<vel_data> vel42_{rc};
  rc_man<mass_data> mass42_{rc};

  rc_man<potential_data> pd42_{rc};

  // Neighbor lists must be initialized after potential initialization.
  // xred, yred, and zred need to be initialized in vdw routines and will be
  // used in nblist setups.
  rc_man<box_data> box42_{rc};
  rc_man<couple_data> cpl42_{rc};
  rc_man<nblist_data> nbl42_{rc};

  rc_man<random_data> rand42_{rc};

  rc_man<md_data> md42_{rc};
}
}
TINKER_NAMESPACE_END
