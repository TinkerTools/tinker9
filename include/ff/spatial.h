#pragma once
#include "precision.h"
#include "tool/genunit.h"

namespace tinker {
struct Spatial
{
   static constexpr int BLOCK = 32;
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

   /// Sorted atoms, containing the x, y, and z coordinates,
   /// and the unsorted atom number.
   struct alignas(16) SortedAtom
   {
      real x, y, z;
      int unsorted;
#if TINKER_REAL_SIZE == 8
      int padding;
#endif
   };
   SortedAtom* sorted; // Length n.
   int* bnum;          // Length n. array[unsorted] == sorted.

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

   ~Spatial();
};
using SpatialUnit = GenericUnit<Spatial, GenericUnitVersion::DISABLE_ON_DEVICE>;

void spatialDataAlloc(SpatialUnit& u, int n, double cutoff, double buffer, //
   const real* x, const real* y, const real* z, int nstype,                //
   int ns1, int (*js1)[2], int ns2, int (*js2)[2],                         //
   int ns3, int (*js3)[2], int ns4, int (*js4)[2]);
void spatialDataInit_cu(SpatialUnit);
void spatialDataUpdateSorted_cu(SpatialUnit);
}
