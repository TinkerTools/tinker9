#pragma once
#include "macro.h"
#include "tool/fc.h"
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
template <class T>
void t_allocate_d1(T** pp, int dim1)
{
   void** p = (void**)pp;
   size_t bytes1 = sizeof(T) * dim1;
   t_allocate_char1(p, bytes1);
}


/**
 * \ingroup bindc
 * \brief Tinker subroutine: version (filename,status).
 * \note This C++ interface does not change the input file string.
 */
std::string t_version(std::string infile, std::string status);
}
