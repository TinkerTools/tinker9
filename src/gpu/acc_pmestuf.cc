#include "acc_e.h"
#include "gpu/acc_fmat.h"
#include "gpu/acc_mathfunc.h"
#include "gpu/decl_box.h"
#include "gpu/decl_mdstate.h"
#include "gpu/decl_pme.h"
#include "gpu/e_mpole.h"

#define TINKER_SRC_GPU_ACC_PMESTUFF_IMPL_
#include "acc_pme_cf.hh"
#include "acc_pme_conv.hh"
#include "acc_pme_grid.hh"
#undef TINKER_SRC_GPU_ACC_PMESTUFF_IMPL_
