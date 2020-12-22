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
 * \brief Tinker function: freeunit ().
 */
int t_freeunit();


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
 * \brief Tinker subroutine: prtxyz
 */
void t_prtxyz(int ixyz);


/**
 * \ingroup bindc
 * \brief Tinker subroutine: prterr
 */
void t_prterr();


//====================================================================//


/**
 * \ingroup bindc
 * \brief Tinker subroutine: optinit
 */
void t_optinit();


/**
 * \ingroup bindc
 * \brief Tinker subroutine: lbfgs (nvar,x0,minimum,grdmin,fgvalue,optsave)
 * \note minimum and optsave are dropped from parameters.
 */
void t_lbfgs(int nvar, double* x0, double grdmin, void* fgvalue);


//====================================================================//


/**
 * \ingroup bindc
 * \brief Tinker subroutine: evcorr1 (mode,elrc,vlrc)
 */
void t_evcorr1(const char* mode, double* elrc, double* vlrc);


//====================================================================//


/**
 * \ingroup bindc
 * \brief Tinker subroutine: extent (rmax)
 */
void t_extent(double& rmax);


//====================================================================//


/**
 * \ingroup bindc
 * \brief Tinker subroutine: lattice
 */
void t_lattice();
TINKER_EXTERN_C_END
