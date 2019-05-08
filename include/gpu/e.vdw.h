#ifndef TINKER_GPU_E_VDW_H_
#define TINKER_GPU_E_VDW_H_

#include "defines.h"

TINKER_NAMESPACE_BEGIN
namespace gpu {
const int evdw_lj = 0x0001;
const int evdw_buck = 0x0002;
const int evdw_mm3hb = 0x0004;
const int evdw_hal = 0x0008;
const int evdw_gauss = 0x0010;
extern int vdwtyp;
extern std::string vdwtyp_str;

const int evdw_decouple = 0;
const int evdw_annihilate = 1;

extern real* vscales;

extern int *jvdw, *njvdw;
extern int* ired;
extern real* kred;
extern real *xred, *yred, *zred;
extern real *radmin, *epsilon;

extern real* ev;
int use_evdw();
void get_evdw_type(int& typ, std::string& typ_str);
real get_evdw();
void e_vdw_data(int op);
}
TINKER_NAMESPACE_END

extern "C" {
void tinker_gpu_evdw_lj0();
void tinker_gpu_evdw_lj1();
void tinker_gpu_evdw_lj3();
void tinker_gpu_evdw_lj4();
void tinker_gpu_evdw_lj5();
void tinker_gpu_evdw_lj6();

void tinker_gpu_evdw_buck0();
void tinker_gpu_evdw_buck1();
void tinker_gpu_evdw_buck3();
void tinker_gpu_evdw_buck4();
void tinker_gpu_evdw_buck5();
void tinker_gpu_evdw_buck6();

void tinker_gpu_evdw_mm3hb0();
void tinker_gpu_evdw_mm3hb1();
void tinker_gpu_evdw_mm3hb3();
void tinker_gpu_evdw_mm3hb4();
void tinker_gpu_evdw_mm3hb5();
void tinker_gpu_evdw_mm3hb6();

void tinker_gpu_evdw_hal0();
void tinker_gpu_evdw_hal1();
void tinker_gpu_evdw_hal3();
void tinker_gpu_evdw_hal4();
void tinker_gpu_evdw_hal5();
void tinker_gpu_evdw_hal6();

void tinker_gpu_evdw_gauss0();
void tinker_gpu_evdw_gauss1();
void tinker_gpu_evdw_gauss3();
void tinker_gpu_evdw_gauss4();
void tinker_gpu_evdw_gauss5();
void tinker_gpu_evdw_gauss6();

void tinker_gpu_evdw0();
void tinker_gpu_evdw1();
void tinker_gpu_evdw3();
void tinker_gpu_evdw4();
void tinker_gpu_evdw5();
void tinker_gpu_evdw6();
}

#endif
