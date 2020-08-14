#pragma once
#include "spatial.h"
#include <map>


namespace tinker {
template <class T>
class Exclusion
{
private:
   static constexpr int MUL = 1000000;
   std::map<int, int> scale_to_enum;
   std::map<int, double> enum_to_scale;

public:
   enum : int
   {
      a = 0,
      b = 1,
      c = 2,
      d = 3,
      e = 4,
      f = 5,
      g = 6,
      h = 7,
   };


   void add_scale_factor(double s)
   {
      if (s < 0 or s > 1)
         return;


      int si = MUL * s;
      auto it = scale_to_enum.find(si);
      if (it == scale_to_enum.end()) {
         // s is a new scale factor
         int sz = scale_to_enum.size();
         scale_to_enum[si] = sz;
         enum_to_scale[sz] = s;
      }
   }


   void clear()
   {
      scale_to_enum.clear();
      enum_to_scale.clear();
   }


   void init()
   {
      add_scale_factor(0.0);
      add_scale_factor(1.0);
   }
};


struct Spatial2
{
   // output


   // internal
   int n, nak;
   int px, py, pz;


   // output


   // internal
   Spatial::SortedAtom* sorted; // n
   int* bnum;                   // n


   int rebuild;
   real cutoff, buffer;
   const real* x;
   const real* y;
   const real* z;
   int* update;              // 2*n
   real *xold, *yold, *zold; // n


   ~Spatial2();
};
using Spatial2Unit = GenericUnit<Spatial2, GenericUnitVersion::EnableOnDevice>;


void spatial2_cut(int& px, int& py, int& pz, int level);
// void* excl_scale ==> real (*excl_scale)[NS]
void spatial2_data_alloc(Spatial2Unit& u, int n, double cutoff, double buffer,
                         const real* x, const real* y, const real* z, int nexcl,
                         int (*excl)[2], void* excl_scale, int NS);
void spatial_data_init_cu(Spatial2Unit);
void spatial_data_update_sorted(Spatial2Unit);


extern Spatial2Unit cspatial_v2_unit;
}
