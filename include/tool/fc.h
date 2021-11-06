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


/**
 * \ingroup bindc
 * \brief Fortran statement: allocated (p).
 */
int t_allocated(void* p);


/**
 * \ingroup bindc
 * \brief Fortran statement: deallocate (p).
 */
void t_deallocate(void* p);


/**
 * \ingroup bindc
 * \brief Fortran statement: allocate (p(dim1)).
 * \code{.f}
 * character, allocatable :: p(:)
 * allocate (p(bytes1))
 * \endcode
 * \note This C interface only allocates bytes (not array elements) in Fortran.
 */
void t_allocate_char1(void** pp, size_t bytes1);


//====================================================================//


/**
 * \ingroup bindc
 * \brief Fortran statement: read (input,10) line, where input is `stdin`,
 * format 10 is format of a character string and line is a character string
 * long enough for format 10.
 */
const char* t_read_stdin_line();


/**
 * \ingroup bindc
 * \brief Tinker subroutine: version (filename,status)
 * \note This C interface does not change the input file string.
 */
const char* t_version(const char* infile, const char* status);


/**
 * \ingroup bindc
 * \brief Tinker subroutine: suffix(filename,extension,status)
 */
void t_suffix(char* filename, const char* extension, const char* status);


/**
 * \ingroup bindc
 * \brief Tinker subroutine: basefile(string)
 */
void t_basefile(char* string);


/**
 * \ingroup bindc
 * \brief Tinker subroutine: prterr
 */
void t_prterr();
TINKER_EXTERN_C_END
