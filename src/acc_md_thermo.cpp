#include "mathfunc.h"
#include "md.h"
#include "random.h"
#include <tinker/detail/bath.hh>
#include <tinker/detail/mdstuf.hh>


TINKER_NAMESPACE_BEGIN
void bussi_thermostat_acc(time_prec dt_prec, T_prec temp_prec)
{
   double dt = dt_prec;
   double temp = temp_prec;


   double tautemp = bath::tautemp;
   double kelvin = bath::kelvin;
   int nfree = mdstuf::nfree;
   double& eta = bath::eta;


   if (temp == 0)
      temp = 0.1;


   double c = std::exp(-dt / tautemp);
   double d = (1 - c) * (kelvin / temp) / nfree;
   double r = normal<double>();
   double s = chi_squared<double>(nfree - 1);
   double scale = c + (s + r * r) * d + 2 * r * std::sqrt(c * d);
   scale = std::sqrt(scale);
   if (r + std::sqrt(c / d) < 0)
      scale = -scale;
   eta *= scale;


   vel_prec sc = scale;
   #pragma acc parallel loop independent async deviceptr(vx,vy,vz)
   for (int i = 0; i < n; ++i) {
      vx[i] *= sc;
      vy[i] *= sc;
      vz[i] *= sc;
   }
}
TINKER_NAMESPACE_END
