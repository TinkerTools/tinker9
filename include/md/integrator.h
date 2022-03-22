#pragma once
#include "ff/rattle.h"
#include "precision.h"
#include <cstdio>
#include <string>
#include <tinker/detail/bath.hh>

namespace tinker {
/// \ingroup mdintg
class IntegratorStaticData
{
protected:
   static int nrespa;
   static bool applyBaro;
   static bool atomic;
   static bool aniso;
   static bool semiiso;
   static double dofP;
   static int anisoArrayLength;
   static const int anisoArray[6][2];
   static constexpr int SemiIso = 2;
   static constexpr int OrthoOrOct = 3;
   static constexpr int Mono = 4;
   static constexpr int Tri = 6;
};

//====================================================================//

/// \ingroup mdintg
enum class PropagatorEnum
{
   Respa,
   Verlet,

   m_LogV,
};

class BasicPropagator;
BasicPropagator* create(PropagatorEnum pe);

/// \ingroup mdintg
/// \brief The interface class of a Verlet or an RESPA-Verlet MD step.
class BasicPropagator : virtual public IntegratorStaticData
{
public:
   /// \brief Logical flag governing saving an MD step.
   /// \param istep  Current number of MD step, started from 1.
   bool ifSave(int istep) const;
   BasicPropagator();
   virtual ~BasicPropagator();

   /// \brief Position update.
   /// \param t  Actual time interval.
   virtual void updatePosition(time_prec t);

   /// \brief The first half-step velocity update.
   /// \param t  Actual time interval, i.e., half-step.
   virtual void updateVelocity1(time_prec t);
   /// \brief The second half-step velocity update.
   /// \param t  Actual time interval, i.e., half-step.
   virtual void updateVelocity2(time_prec t);

   /// \brief Velocity update for the inner RESPA time steps.
   /// \param t  Actual time interval, i.e., time-step/nrespa.
   virtual void updateVelocityR0(time_prec t);
   /// \brief The first half-step velocity update for RESPA.
   /// \param t  Actual time interval, i.e., half-step.
   /// \param nrespa  Number of inner RESPA steps per time-step.
   virtual void updateVelocityR1(time_prec t, int nrespa);
   /// \brief The second half-step velocity update for RESPA.
   /// \param t       Actual time interval, i.e., half-step.
   /// \param nrespa  Number of inner RESPA steps per time-step.
   virtual void updateVelocityR2(time_prec t, int nrespa);

   virtual void rattleSave();
   virtual void rattle(time_prec dt);
   virtual void rattle2(time_prec dt, bool useVirial);
};

typedef BasicPropagator VerletDevice;

class RespaDevice : public BasicPropagator
{
public:
   RespaDevice();
   ~RespaDevice();

   void updateVelocityR1(time_prec t, int nrespa) override;
   void updateVelocityR2(time_prec t, int nrespa) override;
};

class LogVDevice : public BasicPropagator
{
protected:
   typedef BasicPropagator base_t;

   RespaDevice* m_respa;

   void updateVelocityImpl(time_prec t, int idx, int nrespa);

public:
   LogVDevice(bool isNRespa1);
   ~LogVDevice();

   void updatePosition(time_prec t) override;
   void updateVelocity1(time_prec t) override;
   void updateVelocity2(time_prec t) override;

   void updateVelocityR0(time_prec t) override;
   void updateVelocityR1(time_prec t, int nrespa) override;
   void updateVelocityR2(time_prec t, int nrespa) override;
};

//====================================================================//
// generic thermostat

enum class ThermostatEnum
{
   Null,
   Andersen,  // empty
   Berendsen, // empty
   Bussi,
   Nhc,

   m_LeapFrogLP,
   m_LP2022,
   m_Nhc1996,
   m_Nhc2006,
};

class BasicThermostat;
BasicThermostat* create(ThermostatEnum);

class BasicThermostat : virtual public IntegratorStaticData
{
protected:
   void printBasic(FILE*);

public:
   BasicThermostat();
   virtual ~BasicThermostat();
   virtual void printDetail(FILE*);
   virtual void control1(time_prec timeStep);
   virtual void control2(time_prec timeStep, bool save);
};

class BussiThermostat : public BasicThermostat
{
public:
   void printDetail(FILE*) override;
   void control2(double timeStep, bool) override;
};

/// \brief Maximum length of the NH chain.
constexpr int maxnose = bath::maxnose;

/// \ingroup mdpt
/// \brief Applies a velocity correction as needed for the Nose-Hoover Chains
/// at the half time step.
///
/// Literature reference:
///    - <a href="https://doi.org/10.1080/00268979600100761">
///    G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
///    "Explicit Reversible Integrators for Extended Systems Dynamics",
///    Molecular Physics, 87, 1117-1157 (1996).
///    </a>
class NhcDevice : public BasicThermostat
{
protected:
   static constexpr int nhc_nsy = 3;
   int nnose, nhc_nc;
   double g0;
   double vnh[maxnose], qnh[maxnose];
   double (*f_kin)();
   void (*scale_vel)(double);
   std::string name;

   void controlImpl(double timeStep);

public:
   NhcDevice(int nhclen, int nc, double dfree, //
      double (*kin)(), void (*scale)(double), std::string str);
   void printDetail(FILE*);
   void control1(time_prec time_prec) override;
   void control2(time_prec time_prec, bool) override;

   static double kineticAtomic();
   static void scaleVelocityAtomic(double scale);
};

class Nhc06Thermostat : public BasicThermostat
{
protected:
   NhcDevice* m_tpart;
   NhcDevice* m_tbaro;

public:
   Nhc06Thermostat();
   ~Nhc06Thermostat();

   void printDetail(FILE*) override;
   void control1(time_prec dt) override;
   void control2(time_prec dt, bool save) override;

   static double kineticRattleGroup();
   static void scaleVelocityRattleGroup(double scale);

   static double kineticVbar();
   static void scaleVelocityVbar(double scale);

   static double dofVbar();
};

//====================================================================//
// generic barostat

enum class BarostatEnum
{
   Null,
   Berendsen,
   Bussi, // empty
   LP2022,
   MonteCarlo,
   Nhc2006,

   m_LeapFrogLP,
   m_LogVIso,
   m_LogVAniso,
   m_Nhc1996,
};

class BasicBarostat;
BasicBarostat* create(BarostatEnum);

class BasicBarostat : virtual public IntegratorStaticData
{
protected:
   int m_nbaro;
   void printBasic(FILE*);

public:
   BasicBarostat();
   virtual ~BasicBarostat();
   virtual void printDetail(FILE*);
   virtual BarostatEnum getBarostatEnum() const;

   virtual void control1(time_prec timeStep) {}
   virtual void control2(time_prec timeStep) {}
   virtual void control3(time_prec timeStep) {}
   virtual void control4(time_prec timeStep) {}

   virtual bool ifApply(int istep);
   bool ifApply() const;
};

class MonteCarloBarostat : public BasicBarostat
{
public:
   ~MonteCarloBarostat();
   MonteCarloBarostat();
   void printDetail(FILE*) override;
   BarostatEnum getBarostatEnum() const override;
   void control4(time_prec) override;
   bool ifApply(int istep) override;
};

class BerendsenBarostat : public BasicBarostat
{
public:
   void printDetail(FILE*) override;
   BarostatEnum getBarostatEnum() const override;
   void control2(time_prec timeStep) override;
};

class IsoBaroDevice : public BasicBarostat
{
protected:
   double* m_vir;
   double* m_eksum;

   double m_fric;
   double m_rdn;
   bool m_langevin;

   void control_1_2(time_prec dt);

public:
   IsoBaroDevice(double fric);
   BarostatEnum getBarostatEnum() const override;
   void printDetail(FILE*) override;
   void control1(time_prec dt) override;
   void control2(time_prec dt) override;
   void control3(time_prec dt) override;
};

class AnisoBaroDevice : public BasicBarostat
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
   AnisoBaroDevice(double fric);
   BarostatEnum getBarostatEnum() const override;
   void printDetail(FILE*) override;
   void control1(time_prec dt) override;
   void control2(time_prec dt) override;
   void control3(time_prec dt) override;
};

class Nhc06Barostat : public BasicBarostat
{
protected:
   Nhc06Thermostat* m_thermo;
   IsoBaroDevice* m_baro;

public:
   Nhc06Barostat();
   ~Nhc06Barostat();
   void printDetail(FILE*) override;
   BarostatEnum getBarostatEnum() const override;
   void control1(time_prec timeStep) override;
   void control2(time_prec timeStep) override;
   void control3(time_prec timeStep) override;
};

class LP22Barostat : public BasicBarostat
{
protected:
   Nhc06Thermostat* m_thermo;
   BasicBarostat* m_baro;

public:
   LP22Barostat();
   ~LP22Barostat();
   void printDetail(FILE*) override;
   BarostatEnum getBarostatEnum() const override;
   void control1(time_prec timeStep) override;
   void control2(time_prec timeStep) override;
   void control3(time_prec timeStep) override;
};

//====================================================================//
// integrator

enum class IntegratorEnum
{
   Beeman,
   Respa,
   Verlet,

   LeapFrogLP,
   LP2022,
   Nhc1996,
   Nhc2006, // J. Phys. A, 39 5629 (2006), https://doi.org/10.1088/0305-4470/39/19/S18
};

class BasicIntegrator : virtual public IntegratorStaticData
{
protected:
   BasicPropagator* m_prop;
   BasicThermostat* m_thermo;
   BasicBarostat* m_baro;

   int vers1;
   bool save;

   virtual void plan(int istep);
   virtual const char* name() const = 0;
   virtual void kickoff() = 0;

public:
   void printDetail(FILE*);
   BasicIntegrator(PropagatorEnum pe, ThermostatEnum te, BarostatEnum be);
   BasicIntegrator();
   virtual ~BasicIntegrator();
   virtual void dynamic(int istep, time_prec dt);
};

class VerletIntegrator : public BasicIntegrator
{
protected:
   const char* name() const override;
   void kickoff() override;

public:
   VerletIntegrator(ThermostatEnum te, BarostatEnum be);
   static void KickOff();
};

class RespaIntegrator : public BasicIntegrator
{
protected:
   const char* name() const override;
   void kickoff() override;

public:
   RespaIntegrator(ThermostatEnum te, BarostatEnum be);
   static void KickOff();
};

class Nhc96Integrator : public BasicIntegrator
{
protected:
   const char* name() const override;
   void kickoff() override;

public:
   Nhc96Integrator();
   void dynamic(int, time_prec) override;
};

class Nhc06Integrator : public BasicIntegrator
{
protected:
   bool m_isNRespa1;
   const char* name() const override;
   void kickoff() override;

public:
   Nhc06Integrator(bool isNRespa1);
};

class LP22Integrator : public BasicIntegrator
{
protected:
   bool m_isNRespa1;
   const char* name() const override;
   void kickoff() override;

public:
   LP22Integrator(bool isNRespa1);
};

class LeapFrogLPIntegrator : public BasicIntegrator
{
protected:
   const char* name() const override;
   void kickoff() override;

public:
   LeapFrogLPIntegrator();
   ~LeapFrogLPIntegrator();
   void dynamic(int, time_prec) override;
};
}
