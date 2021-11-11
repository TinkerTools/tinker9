#pragma once
#include <string.h>


#ifdef __cplusplus
#   define TINKER_EXTERN_C_BEGIN                                               \
      extern "C"                                                               \
      {
#else
#   define TINKER_EXTERN_C_BEGIN
#endif
#ifdef __cplusplus
#   define TINKER_EXTERN_C_END }
#else
#   define TINKER_EXTERN_C_END
#endif


TINKER_EXTERN_C_BEGIN
/**
 * \ingroup bindc
 * \brief Fortran statement: rewind (unit=iunit).
 */
void t_rewind(int unit);


/**
 * \ingroup bindc
 * \brief Fortran statement: close (unit=iunit).
 */
void t_close(int unit);


/**
 * \ingroup bindc
 * \brief Fortran statement: open (unit=iunit,file=ifile,status=istatus).
 */
void t_open(int unit, const char* file, const char* status);


//====================================================================//


//====================================================================//


/**
 * \ingroup bindc
 * \brief Tinker subroutine: prterr
 */
void t_prterr();
TINKER_EXTERN_C_END

#include "macro.h"
#include <string>


namespace tinker {
/**
 * \ingroup bindc
 * \brief Fortran statement: open (unit=iunit,file=ifile,status=istatus).
 */
void t_open(int unit, std::string file, std::string status);


/**
 * \ingroup bindc
 * \brief Fortran statement: allocate (p(dim1)).
 * \note This C++ interface allocates array elements in Fortran.
 */

}


#include <tinker/routines.h>

// version
std::string tinker_f_version(std::string infile, std::string status);

// memory
int tinker_f_allocated(void* p);

void tinker_f_deallocate(void* p);

void tinker_f_allocate_byte(void** pp, size_t bytes);

template <class T>
void tinker_f_allocate_dim(T** pp, int dim)
{
   void** p = (void**)pp;
   size_t bytes = sizeof(T) * dim;
   tinker_f_allocate_byte(p, bytes);
}

// read stdin
std::string tinker_f_read_stdin_line();
