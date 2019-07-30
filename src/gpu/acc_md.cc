#include "mod_mdstate.h"
#include "util_io.h"
#include "util_mdstate.h"
#include <ext/tinker/tinker_mod.h>
#include <ext/tinker/tinker_rt.h>

TINKER_NAMESPACE_BEGIN
void kinetic_acc_impl_(real& temp) {
  const real ekcal_inv = 1.0 / units::ekcal;
  real exx = 0;
  real eyy = 0;
  real ezz = 0;
  real exy = 0;
  real eyz = 0;
  real ezx = 0;
  #pragma acc parallel loop independent deviceptr(mass,vx,vy,vz)\
          reduction(+:exx,eyy,ezz,exy,eyz,ezx)
  for (int i = 0; i < n; ++i) {
    real term = 0.5f * mass[i] * ekcal_inv;
    exx += term * vx[i] * vx[i];
    eyy += term * vy[i] * vy[i];
    ezz += term * vz[i] * vz[i];
    exy += term * vx[i] * vy[i];
    eyz += term * vy[i] * vz[i];
    ezx += term * vz[i] * vx[i];
  }
  ekin[0][0] = exx;
  ekin[0][1] = exy;
  ekin[0][2] = ezx;
  ekin[1][0] = exy;
  ekin[1][1] = eyy;
  ekin[1][2] = eyz;
  ekin[2][0] = ezx;
  ekin[2][1] = eyz;
  ekin[2][2] = ezz;
  eksum = exx + eyy + ezz;
  temp = 2 * eksum / (mdstuf::nfree * units::gasconst);

  // TODO: RIGIDBODY
  // TODO: if (isobaric .and. barostat.eq.'BUSSI')
}

void mdrest_acc_impl_(int istep) {
  if (!mdstuf::dorest)
    return;
  if ((istep % mdstuf::irest) != 0)
    return;

  // TODO: special treatment for 'RIGIDBODY'

  const real ekcal = units::ekcal;

  // zero out the total mass and overall linear velocity

  real totmass = 0;
  real vtot1 = 0;
  real vtot2 = 0;
  real vtot3 = 0;

  // compute linear velocity of the system center of mass

  #pragma acc parallel loop independent deviceptr(mass,vx,vy,vz)\
              reduction(+:totmass,vtot1,vtot2,vtot3)
  for (int i = 0; i < n; ++i) {
    real weigh = mass[i];
    totmass += weigh;
    vtot1 += vx[i] * weigh;
    vtot2 += vy[i] * weigh;
    vtot3 += vz[i] * weigh;
  }

  vtot1 /= totmass;
  vtot2 /= totmass;
  vtot3 /= totmass;

  // compute translational kinetic energy of overall system

  real etrans = vtot1 * vtot1 + vtot2 * vtot2 + vtot3 * vtot3;
  etrans *= 0.5f * totmass / ekcal;

  real erot, xtot, ytot, ztot, vang[3];
  if (!bound::use_bounds) {

    // find the center of mass coordinates of the overall system
    // compute the angular momentum of the overall system

    xtot = 0;
    ytot = 0;
    ztot = 0;

    real mang1 = 0;
    real mang2 = 0;
    real mang3 = 0;

    #pragma acc parallel loop independent deviceptr(mass,x,y,z,vx,vy,vz)\
                reduction(+:xtot,ytot,ztot,mang1,mang2,mang3)
    for (int i = 0; i < n; ++i) {
      real weigh = mass[i];
      xtot += x[i] * weigh;
      ytot += y[i] * weigh;
      ztot += z[i] * weigh;
      mang1 += (y[i] * vz[i] - z[i] * vy[i]) * weigh;
      mang2 += (z[i] * vx[i] - x[i] * vz[i]) * weigh;
      mang3 += (x[i] * vy[i] - y[i] * vx[i]) * weigh;
    }
    xtot /= totmass;
    ytot /= totmass;
    ztot /= totmass;
    mang1 -= (ytot * vtot3 - ztot * vtot2) * totmass;
    mang2 -= (ztot * vtot1 - xtot * vtot3) * totmass;
    mang3 -= (xtot * vtot2 - ytot * vtot1) * totmass;

    // calculate the moment of inertia tensor

    real xx = 0;
    real xy = 0;
    real xz = 0;
    real yy = 0;
    real yz = 0;
    real zz = 0;

    #pragma acc parallel loop independent deviceptr(mass,x,y,z)\
                reduction(+:xx,xy,xz,yy,yz,zz)
    for (int i = 0; i < n; ++i) {
      real weigh = mass[i];
      real xdel = x[i] - xtot;
      real ydel = y[i] - ytot;
      real zdel = z[i] - ztot;
      xx += xdel * xdel * weigh;
      xy += xdel * ydel * weigh;
      xz += xdel * zdel * weigh;
      yy += ydel * ydel * weigh;
      yz += ydel * zdel * weigh;
      zz += zdel * zdel * weigh;
    }

    double tensor[3][3];
    double eps = (n <= 2 ? 0.000001 : 0);
    tensor[0][0] = yy + zz + eps;
    tensor[0][1] = -xy;
    tensor[0][2] = -xz;
    tensor[1][0] = -xy;
    tensor[1][1] = xx + zz + eps;
    tensor[1][2] = -yz;
    tensor[2][0] = -xz;
    tensor[2][1] = -yz;
    tensor[2][2] = xx + yy + eps;

    int ndim = 3;
    TINKER_RT(invert)(&ndim, &tensor[0][0]);

    // compute angular velocity and rotational kinetic energy

    for (int i = 0; i < 3; ++i) {
      vang[i] =
          tensor[0][i] * mang1 + tensor[1][i] * mang2 + tensor[2][i] * mang3;
    }
    erot = vang[0] * mang1 + vang[1] * mang2 + vang[2] * mang3;
    erot *= (0.5f / ekcal);
  }

  // eliminate any translation of the overall system

  #pragma acc parallel loop independent deviceptr(vx,vy,vz)
  for (int i = 0; i < n; ++i) {
    vx[i] -= vtot1;
    vy[i] -= vtot2;
    vz[i] -= vtot3;
  }

  // print the translational velocity of the overall system

  if (inform::debug) {
    print(stdout,
          " System Linear Velocity :  {:12.2f}{:12.2f}{:12.2f}\n Translational "
          "Kinetic Energy :{:10s}{:12.4f} Kcal/mole\n",
          vtot1, vtot2, vtot3, "", etrans);
  }

  // eliminate any rotation about the system center of mass

  if (!bound::use_bounds) {
    #pragma acc parallel loop independent deviceptr(x,y,z,vx,vy,vz)
    for (int i = 0; i < n; ++i) {
      real xdel = x[i] - xtot;
      real ydel = y[i] - ytot;
      real zdel = z[i] - ztot;
      vx[i] -= vang[1] * zdel + vang[2] * ydel;
      vy[i] -= vang[2] * xdel + vang[0] * zdel;
      vz[i] -= vang[0] * ydel + vang[1] * xdel;
    }
  }

  // print the angular velocity of the overall system

  if (inform::debug) {
    print(stdout,
          " System Angular Velocity : {:12.2f}{:12.2f}{:12.2f}\n Rotational "
          "Kinetic Energy :{:13s}{:12.4f} Kcal/mole\n",
          vang[0], vang[1], vang[2], "", erot);
  }
}

//======================================================================
// propagate

void propagate_xyz_acc_impl_(real dt) {
  #pragma acc parallel loop independent deviceptr(x,y,z,vx,vy,vz)
  for (int i = 0; i < n; ++i) {
    x[i] += dt * vx[i];
    y[i] += dt * vy[i];
    z[i] += dt * vz[i];
  }
}

void propagate_velocity_acc_impl_(real dt) {
  const real ekcal = units::ekcal;
  #pragma acc parallel loop independent deviceptr(vx,vy,vz,gx,gy,gz,massinv)
  for (int i = 0; i < n; ++i) {
    real coef = -ekcal * massinv[i] * dt;
    vx[i] += coef * gx[i];
    vy[i] += coef * gy[i];
    vz[i] += coef * gz[i];
  }
}
TINKER_NAMESPACE_END
