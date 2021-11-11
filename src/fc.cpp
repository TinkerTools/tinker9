#include "tool/fc.h"
#include "box.h"
#include "mdpq.h"
#include "tool/darray.h"
#include <tinker/detail/atoms.hh>
#include <tinker/routines.h>


namespace {
char out[2048];
static_assert(2048 >= tinker::MAX_NCHAR + 5, "");
}


TINKER_EXTERN_C_BEGIN
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


extern "C" void fc_allocated_(void** p, int* ans);
int t_allocated(void* p)
{
   int ans;
   fc_allocated_(&p, &ans);
   return ans;
}


extern "C" void fc_deallocate_(void**);
void t_deallocate(void* p)
{
   fc_deallocate_(&p);
}


extern "C" void fc_allocate_char1_(void**, size_t);
void t_allocate_char1(void** pp, size_t bytes1)
{
   fc_allocate_char1_(pp, bytes1);
}


//====================================================================//


extern "C" void fc_read_stdin_line(char* outfile);
const char* t_read_stdin_line()
{
   fc_read_stdin_line(out);
   return out;
}


extern "C" void fc_version(char* outfile, const char* file, int flen,
                           const char* status, int slen);
const char* t_version(const char* file, const char* status)
{
   int flen = strlen(file);
   int slen = strlen(status);
   fc_version(out, file, flen, status, slen);
   return out;
}


extern "C" void fc_suffix(char* file, const char* ext, const char* status,
                          int slen);
void t_suffix(char* filename, const char* extension, const char* status)
{
   int slen = strlen(status);
   fc_suffix(filename, extension, status, slen);
}


void t_prterr()
{
   using namespace tinker;
   Box p;
   get_default_box(p);
   set_tinker_box_module(p);
   darray::copyout(g::q0, n, atoms::x, xpos);
   darray::copyout(g::q0, n, atoms::y, ypos);
   darray::copyout(g::q0, n, atoms::z, zpos);
   wait_for(g::q0);
   tinker_f_prterr();
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


std::string t_read_stdin_line()
{
   fc_read_stdin_line(out);
   return std::string(out);
}


std::string t_version(std::string file, std::string status)
{
   int flen = file.size();
   int slen = status.size();
   fc_version(out, file.c_str(), flen, status.c_str(), slen);
   return std::string(out);
}
}

extern "C"
{}
