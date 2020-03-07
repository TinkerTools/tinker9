#include "evdw.h"


TINKER_NAMESPACE_BEGIN
template <class Ver>
void elj_acc1()
{
   constexpr bool do_e = Ver::e;
   constexpr bool do_a = Ver::a;
   constexpr bool do_g = Ver::g;
   constexpr bool do_v = Ver::v;
}


void elj_acc(int vers) {}


void ebuck_acc(int) {}
void emm3hb_acc(int) {}
void egauss_acc(int) {}
TINKER_NAMESPACE_END
