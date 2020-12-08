Multipole Moment
================

.. include:: ../replace.rst

Definitions and Units
---------------------

The electrostatic potential at *r* due to the charge distribution nearby is

.. math::
   :label: pot1

   \phi(r) = \frac{1}{4\pi\epsilon_0}
   \int ds\frac{\rho(s)}{|r-s|},\ ds=\mathrm{d}x \mathrm{d}y \mathrm{d}z.

Tinker uses a variable **electric** (in **chgpot** module) to represent the
the factor :math:`1/(4\pi\epsilon_0)`.
Its default magnitude is 332.063713, which is a constant defined by variable
**coulomb** (in **units** module), and its units are kcal/mol |ang|/|e2|.
The default value can be modified by the *ELECTIRIC* keyword.

.. note::

   Should the value of **coulomb** documented here be out-dated and become
   inconsistent with our code, please send us a pull request.

Expanding :math:`1/|r-s|` in Taylor series, :math:`4\pi\epsilon_0\phi(r)`
can be rewritten as

.. math::
   :label: pot2

   \left[\int ds\rho(s)\right]\frac{1}{r}
   -\sum_i\left[\int ds\rho(s)s_i\right]\nabla_i\frac{1}{r}
   +\sum_{ij}\left[\frac{1}{2}\int ds\rho(s)s_i s_j\right]\nabla_i\nabla_j\frac{1}{r}
   - \cdots,

where three pairs of square brackets give a set of definitions of monopole
(charge, *C*), dipole (*D*), and quadrupole moments (*Q\**), respectively.
The units of the multipole moments used in Tinker parameter files and internal
calculation are different.

==========  ===============  ==============
Multipole   Parameter Units  Internal Units
==========  ===============  ==============
Charge      e                e
Dipole      e Bohr           e |ang|
Quadrupole  e |bohr2|        e |ang2|
==========  ===============  ==============

In addition to different units, the quadrupole moments in Tinker parameter
files use what is traditionally called *traceless quadrupole* :math:`\Theta`
that has a different definition than *Q\**.
The third term in :eq:`pot2` can be rewritten as

.. math::

   \sum_{ij}\left[\frac{1}{2}\int ds\rho(s)(3s_i s_j - s^2\delta_{ij})\right]
   \frac{r_i r_j}{r^5},

hence the traceless quadrupole can be defined as

.. math::

   \Theta_{ij} = \frac{1}{2}\int ds\rho(s)(3s_i s_j - s^2\delta_{ij}).

It is easy to confirm that :math:`\sum_k^{x,y,z}(3 s_k s_k - s^2)=0`, thus,

.. math::

   \Theta_{ij} = 3Q_{ij}^* - \delta_{ij}\sum_k^{x,y,z}Q_{kk}^*.

Internally, Tinker scales :math:`\Theta` by 1/3

.. math::

   Q = \Theta/3,

so that the energy expression is the same as if we were using *Q\**.

Energy Torque  Gradient
-----------------------

Potential energy

.. math::
   :label: pot3

   U = \frac{1}{4\pi\epsilon_0}\int ds\rho(s)\phi(s).

Potential energy with discretized charge distribution in :eq:`pot3`

.. math::
   :label: pot4

   U(r) = \phi(r) C(r) + \phi'(r) D(r) + \phi''(r) Q(r) + \cdots.

Distance

.. math::

   (r_x,r_y,r_z)=\boldsymbol{r}=r_2-r_1.

Pairwise (atoms 1 and 2) quadrupole energy

.. math::

   U_{12} = M_1^T T_{12} M_2.

Multipoles

.. math::

   M_1 = \begin{pmatrix}
   C_1 \\
   D_1 \\
   Q_1 \end{pmatrix},\ M_2 = \begin{pmatrix}
   C_2 \\
   D_2 \\
   Q_2 \end{pmatrix}.

T matrix

.. math::

   T_{12} = \begin{pmatrix}
   1          & \nabla_2           & \nabla_2^2           \\
   \nabla_1   & \nabla_1\nabla_2   & \nabla_1\nabla_2^2   \\
   \nabla_1^2 & \nabla_1^2\nabla_2 & \nabla_1^2\nabla_2^2
   \end{pmatrix}\frac{1}{r}.

The upper left 4\ |x|\ 4 elements of :math:`T_{12}`

.. math::

   T_{12}^{4 \times 4} = \begin{pmatrix}
   1/r     & -r_x/r^3           & -r_y/r^3           & -r_z/r^3     \\
   r_x/r^3 & -3r_x^2/r^5 +1/r^3 & -3r_xr_y/r^5       & -3r_xr_z/r^5 \\
   r_y/r^3 & -3r_xr_y/r^5       & -3r_y^2/r^5 +1/r^3 & -3r_yr_z/r^5 \\
   r_z/r^3 & -3r_xr_z/r^5       & -3r_yr_z/r^5       & -3r_z^2/r^5 +1/r^3
   \end{pmatrix}.

In the EWALD summation, :math:`1/r^k` terms will have different forms (*Bn*).
Neverthelss, they are still connected through derivatives.

=================================  =====================================
Non-EWALD                          EWALD
=================================  =====================================
:math:`1/r`                        :math:`B_0=\mathrm{erfc}(\alpha r)/r`
:math:`1/r^3`                      :math:`B_1`
:math:`3/r^5`                      :math:`B_2`
:math:`15/r^7`                     :math:`B_3`
:math:`105/r^9`                    :math:`B_4`
:math:`945/r^{11}`                 :math:`B_5`
=================================  =====================================

The *Bn* terms are related to the (complementary) Boys functions and
(complementary) error functions. For :math:`x>0` and :math:`n\ge 0`,

.. math::
   \frac{\mathrm{erf}(x)}{x} &= \frac{2}{\sqrt{\pi}}F_0(x^2),    \\
   \frac{\mathrm{erfc}(x)}{x} &= \frac{2}{\sqrt{\pi}}F_0^C(x^2), \\
   F_n(x) &= \int_0^1 \exp(-xt^2)t^{2n} dt,                      \\
   F_n^C(x) &= \int_1^\infty \exp(-xt^2)t^{2n} dt.

The Boys functions can be generated through upward and downward recursions

.. math::

   F_n(x) = \frac{2xF_{n+1}(x) + \exp(-x)}{2n+1}, \\
   F_n^C(x) = \frac{2xF_{n+1}^C(x) - \exp(-x)}{2n+1}.

Energy, torque, and force

=========  ================  ======================  ====================
Terms      Energy            Torque                  Force
=========  ================  ======================  ====================
C          :math:`\phi C`    N/A                     :math:`\phi' C`
D          :math:`\phi' D`   :math:`\phi'\times D`   :math:`\phi'' D`
Q          :math:`\phi'' Q`  :math:`\phi''\times Q`  :math:`\phi''' Q`
=========  ================  ======================  ====================

.. math::

   \tau(D) &= \phi'\times D = D \times E, \\
   \tau_{i}(Q) &= -2 \sum_{jk}\sum_{l}^{xyz} \epsilon_{ijk}Q_{jl}\phi''_{kl},

where :math:`\epsilon_{ijk}` is the Levi-Civita symbol.

Reference: :cite:`mpole-pme`.
