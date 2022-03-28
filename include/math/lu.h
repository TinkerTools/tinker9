#pragma once

namespace tinker {
/// \ingroup math
/// \brief This subroutine uses the LU decomposition method to solve the linear
/// system Ax = b, returning x in b. A is an n by n real symmetric matrix with
/// its upper triangle (including the diagonal) stored by rows.
///
/// Literature reference:
///    - <a href="http://numerical.recipes">
///    W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery,
///    Numerical Recipes (C++), 3rd Ed., Section 2.3,
///    Cambridge University Press (2007).
///    </a>
template <int n, class R>
void symlusolve(const R* aUpRowMajor, R* b);
}
