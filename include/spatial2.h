#pragma once
#include "spatial.h"


namespace tinker {
/**
 * \ingroup spatial2
 * \page spatial2
 *
 * %Spatial Decomposition Version 2.
 *
 * ### A. Concepts
 *    - Number of atoms (`n`).
 *    - Number of atom blocks (`nak`), (n+32-1)/32.
 *    - Number of atom block pairs (`nakp`), nak*(nak+1)/2.
 *    - SortedAtomBlockPairs (`akp`).
 *       - For any two atom blocks `(x,y)`, where
 *       `0 <= y <= x < nak`, we can map them to `(row,column)` of a triangle
 *       to obtain a unique number `f(x,y)`.
 *    ```
 *    0
 *    1 2
 *    3 4 5
 *    6 7 8 9
 *    ...
 *    ```
 *       - It is obvious that \f$ x(x+1)/2 \le f(x,y) < (x+1)(x+2)/2 \f$.
 *       `f` is represented by a signed 32-bit integer in the program, thus its
 *       upper limit is \f$ 2^{32}-1 \f$, and the upper limit of `x` is 65535.
 *       The number of atoms should then not exceed 32*(65535+1), i.e. 2097152.
 *       Using the left side of the inequality, we have
 *       \f$ x = \lfloor( \sqrt{8f+1} - 1)/2 \rfloor \f$, where double precision
 *       floating-point arithmetic will be sufficient.
 *    - AtomBlockPairFlags (`akpf`): 32-bit integer array of length `nakpk`,
 *    which equals (nakp+32-1)/32. The i-th bit is set if the i-th
 *    SortedAtomBlockPair is recorded in the neighbor list.
 *    - Boxes: The periodic boundary box ("the big box") is partitioned
 *    into smaller boxes.
 *       - The ranges of the fractional coordinates are all `[-0.5, 0.5)`.
 *       Along x, y, and z axes, each direction is equally split into
 *       \f$ 2^p \f$, \f$ 2^p \f$, and \f$ 2^p \f$ parts, respectively.
 *       - Every box can be accessed by 3 integers `(ix,iy,iz)`, all of which
 *       have range \f$ [0,2^p) \f$. These boxes are placed on a 3-D Hilbert
 *       curve so that we can map `(ix,iy,iz)` to a unique integer value
 *       `BoxID`, which has the range \f$ [0,2^{3p}) \f$. Since `BoxID` is
 *       a signed integer in the program, \f$ p \le 10 \f$.
 *       - <a href="https://doi.org/10.1063/1.1751381">
 *       John Skilling, "Programming the Hilbert curve",
 *       AIP Conf. Proc., 707, (2004).
 *       </a>
 *       - <a href="https://stackoverflow.com/a/10384110">
 *       Paul Chernoch (stackoverflow question 499166, answer 10384110)
 *       </a>
 *       corrected the typo in `TransposetoAxes()`, which should read
 *    ```
 *    for(i=n-1; i>0; i--) X[i] ^= X[i-1];
 *    ```
 *
 * ### B. Step 1
 *    1. Based on the fractional coordinates, assign every atom to a box,
 *    save `(BoxID,AtomNum)` in a temporary array (`b2num`) of size `(2,n)`.
 *    2. Zero out `akpf`.
 *    3. Sort `b2num` by `BoxID`.
 *
 * ### C. Step 2
 *    1. For every sorted atom, save `SortedAtomNum` in `bnum[AtomNum]`.
 *    2. Save sorted coordinates.
 *    3. Zero out `b2num`.
 *    4. For each atom block, compute mid point, radius, half size, and the
 *    "local flag".
 *
 *
 * ### D. Step 3
 *    1. For every atom block, set bit in `akpf` for block pair `(i,i)`.
 */


/**
 * \ingroup spatial2
 */
struct Spatial2
{
   static constexpr int LSTCAP = 48;
   // output
   int nakpl;  // Length of iakpl.
   int niak;   // Length of iak. Not greater than LSTCAP*nak.
   int* iakpl; // List of recorded block pairs. Length nakpl.
   int* iak;   // List of blocks in "block-atoms".
   int* lst;   // List of atoms in "block-atoms".


   // internal
   int n;     // Number of atoms.
   int nak;   // Number of blocks.
   int nakp;  // Number of block pairs. (nak+1)*nak/2.
   int nakpk; // Length of 32-bit integer array to store atom block pairs.
   int px, py, pz;
   int cap_nakpl; // Capacity of iakpl. Initial value (32+8*nak).


   // internal
   int* iakpl_rev; // Length nakp. array[pair] == location in iakpl.
   int* akpf;      // Length nakpk. Block pair bit flags.


   Spatial::SortedAtom* sorted; // Length n.
   int* bnum;                   // Length n. array[unsorted] == sorted.


   struct alignas(16) Center
   {
      real x, y, z, w;
   };
   Center* akc;  // Length nak. Block center and the "local flag".
   Center* half; // Length nak. Half box size and radius.


   int fresh;
   real cutoff, buffer;
   const real* x;
   const real* y;
   const real* z;
   int* update;              // Length max(2*n,128).
   real *xold, *yold, *zold; // Length n.


   struct ScaleInfo
   {
      int (*js)[2];       // Length ns. Atom pairs.
      unsigned int* bit0; // Length 32*cap_nakpl.
      int ns;

      void init();
      void set(int nns, int (*jjs)[2]);
   };
   int nstype; // number of ScaleInfo objects in-use
   ScaleInfo si1;
   ScaleInfo si2;
   ScaleInfo si3;
   ScaleInfo si4;


   ~Spatial2();
};
using Spatial2Unit = GenericUnit<Spatial2, GenericUnitVersion::DisableOnDevice>;


void spatial2_data_alloc(Spatial2Unit& u, int n, double cutoff, double buffer,
                         const real* x, const real* y, const real* z,
                         int nstype,                                     //
                         int ns1, int (*js1)[2], int ns2, int (*js2)[2], //
                         int ns3, int (*js3)[2], int ns4, int (*js4)[2]);
void spatial2_cut(int& px, int& py, int& pz, int level);
void spatial_data_init_cu(Spatial2Unit);
void spatial_data_update_sorted(Spatial2Unit);
}
