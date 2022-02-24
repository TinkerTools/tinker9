#pragma once
#include "itgbBasic.h"

namespace tinker {
class LogVAnisoBarostat : public BasicBarostat
{
protected:
   double* m_vir;
   double* m_eksum;
   double (*m_ekin)[3];

   double m_fric;
   double m_rdn[3][3];
   bool m_langevin;

   void control_1_2(time_prec dt);

public:
   LogVAnisoBarostat(double fric);
   BarostatEnum getBarostatEnum() const override;
   void printDetail(FILE*) override;
   void control1(time_prec dt) override;
   void control2(time_prec dt) override;
   void control3(time_prec dt) override;
};
}
