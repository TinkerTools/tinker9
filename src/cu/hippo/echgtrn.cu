#include "ff/elec.h"
#include "ff/hippomod.h"
#include "ff/image.h"
#include "ff/spatial.h"
#include "ff/switch.h"
#include "seq/add.h"
#include "seq/launch.h"
#include "seq/pair_chgtrn.h"
#include "seq/triangle.h"
#include <cassert>

namespace tinker {
#include "echgtrn_cu1.cc"

template <class Ver>
void echgtrn_cu2()
{
   const auto& st = *mspatial_v2_unit;
   real cut = switchCut(Switch::CHGTRN);
   real off = switchOff(Switch::CHGTRN);
   real f = electric / dielec;

   assert(ctrntyp == Chgtrn::SEPARATE);
   int ngrid = gpuGridSize(BLOCK_DIM);
   echgtrn_cu1<Ver><<<ngrid, BLOCK_DIM, 0, g::s0>>>(st.n, TINKER_IMAGE_ARGS, nct, ect, vir_ect, dectx, decty, dectz,
      cut, off, st.si1.bit0, nmdwexclude, mdwexclude, mdwexclude_scale, st.x, st.y, st.z, st.sorted, st.nakpl, st.iakpl,
      st.niak, st.iak, st.lst, chgct, dmpct, f);
}

void echgtrn_cu(int vers)
{
   if (vers == calc::v0)
      echgtrn_cu2<calc::V0>();
   else if (vers == calc::v1)
      echgtrn_cu2<calc::V1>();
   else if (vers == calc::v3)
      echgtrn_cu2<calc::V3>();
   else if (vers == calc::v4)
      echgtrn_cu2<calc::V4>();
   else if (vers == calc::v5)
      echgtrn_cu2<calc::V5>();
   else if (vers == calc::v6)
      echgtrn_cu2<calc::V6>();
}
}
