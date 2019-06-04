#include "gpu/decl_mdstate.h"
#include "gpu/e_polar.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
template <int USE>
void epolar_coulomb_tmpl() {}
}
TINKER_NAMESPACE_END

extern "C" {
m_tinker_using_namespace;
void tinker_gpu_epolar_coulomb0() { gpu::epolar_coulomb_tmpl<gpu::v0>(); }
void tinker_gpu_epolar_coulomb1() { gpu::epolar_coulomb_tmpl<gpu::v1>(); }
void tinker_gpu_epolar_coulomb3() { gpu::epolar_coulomb_tmpl<gpu::v3>(); }
void tinker_gpu_epolar_coulomb4() { gpu::epolar_coulomb_tmpl<gpu::v4>(); }
void tinker_gpu_epolar_coulomb5() { gpu::epolar_coulomb_tmpl<gpu::v5>(); }
void tinker_gpu_epolar_coulomb6() { gpu::epolar_coulomb_tmpl<gpu::v6>(); }
}
