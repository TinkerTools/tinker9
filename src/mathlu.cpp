namespace tinker {
template <int n, class R>
extern void symlusolve_acc(const R* aUpRowMajor, R* b);
}

namespace tinker {
template <int n, class R>
void symlusolve(const R* aUpRowMajor, R* b)
{
   symlusolve_acc<n, R>(aUpRowMajor, b);
}

template void symlusolve<3, float>(const float*, float*);
template void symlusolve<6, float>(const float*, float*);
template void symlusolve<3, double>(const double*, double*);
template void symlusolve<6, double>(const double*, double*);
}
