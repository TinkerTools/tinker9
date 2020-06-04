#include "tool/fc.h"
#include "tool/fcxx.h"


namespace {
constexpr int MAX_NCHAR = 240;
char out[2048];
static_assert(2048 >= MAX_NCHAR + 5, "");
}


TINKER_EXTERN_C_BEGIN
int TINKER_RT(freeunit)();
int t_freeunit()
{
   return TINKER_RT(freeunit)();
}


void fc_rewind(int);
void t_rewind(int unit)
{
   fc_rewind(unit);
}


void fc_close(int);
void t_close(int unit)
{
   fc_close(unit);
}


void fc_open(int u, const char* f, int fl, const char* st, int sl);
void t_open(int unit, const char* file, const char* status)
{
   int fl = strlen(file);
   int sl = strlen(status);
   fc_open(unit, file, fl, status, sl);
}


//====================================================================//


extern "C" void fc_allocated(void** p, int* ans);
int t_allocated(void* p)
{
   int ans;
   fc_allocated(&p, &ans);
   return ans;
}


extern "C" void fc_deallocate(void**);
void t_deallocate(void* p)
{
   fc_deallocate(&p);
}


extern "C" void fc_allocate_char1(void**, size_t);
void t_allocate_char1(void** pp, size_t bytes1)
{
   fc_allocate_char1(pp, bytes1);
}


//====================================================================//


extern "C" void fc_version(char* outfile, const char* file, int flen,
                           const char* status, int slen);
const char* t_version(const char* file, const char* status)
{
   int flen = strlen(file);
   int slen = strlen(status);
   fc_version(out, file, flen, status, slen);
   return out;
}


extern "C" void TINKER_RT(prtxyz)(int*);
void t_prtxyz(int ixyz)
{
   TINKER_RT(prtxyz)(&ixyz);
}


//====================================================================//


extern "C" void TINKER_RT(optinit)();
void t_optinit()
{
   TINKER_RT(optinit)();
}


extern "C" void TINKER_RT(optsave)(int* ncycle, double* f, double* xx);
extern "C" void TINKER_RT(lbfgs)(int* nvar, double* x0, double* minimum,
                                 double* grdmin, void* fgvalue, void* optsave);
void t_lbfgs(int nvar, double* x0, double grdmin, void* fgvalue)
{
   double minimum;
   void* opts = (void*)TINKER_RT(optsave);
   TINKER_RT(lbfgs)(&nvar, x0, &minimum, &grdmin, fgvalue, opts);
}


//====================================================================//


void fc_evcorr1(const char* mode, int mlen, double* elrc, double* vlrc);
void t_evcorr1(const char* mode, double* elrc, double* vlrc)
{
   int mlen = strlen(mode);
   fc_evcorr1(mode, mlen, elrc, vlrc);
}


//====================================================================//


void TINKER_RT(extent)(double*);
void t_extent(double& rmax)
{
   TINKER_RT(extent)(&rmax);
}


//====================================================================//


void TINKER_RT(lattice)();
void t_lattice()
{
   TINKER_RT(lattice)();
}
TINKER_EXTERN_C_END


//====================================================================//


namespace tinker {
void t_open(int unit, std::string file, std::string status)
{
   int fl = file.size();
   int sl = status.size();
   fc_open(unit, file.c_str(), fl, status.c_str(), sl);
}


std::string t_version(std::string file, std::string status)
{
   int flen = file.size();
   int slen = status.size();
   fc_version(out, file.c_str(), flen, status.c_str(), slen);
   return std::string(out);
}
}
