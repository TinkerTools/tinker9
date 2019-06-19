#include "gpu/decl_box.h"
#include "gpu/decl_dataop.h"
#include "rc_cudart.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
box_t* box;

void box_data(int op) {
  if (op & op_dealloc) {
    check_cudart(cudaFree(box));
  }

  // TODO
  if (op & op_alloc) {
  }
  if (op & op_copyin) {
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

    size_t size = sizeof(box_t);
    check_cudart(cudaMalloc(&box, size));

    copyin_data(&box->xbox, &boxes::xbox, 1);
    copyin_data(&box->ybox, &boxes::ybox, 1);
    copyin_data(&box->zbox, &boxes::zbox, 1);
    copyin_data(&box->alpha, &boxes::alpha, 1);
    copyin_data(&box->beta, &boxes::beta, 1);
    copyin_data(&box->gamma, &boxes::gamma, 1);
    copyin_data(&box->lvec[0][0], &boxes::lvec[0][0], 9);
    copyin_data(&box->recip[0][0], &boxes::recip[0][0], 9);
    copyin_data(&box->volbox, &boxes::volbox, 1);
    copyin_data(&box->shape, &shape, 1);
  }
}
}
TINKER_NAMESPACE_END
