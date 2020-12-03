Multipole Moment
================

.. include:: ../replace.rst

The electrostatic potential at *r* due to the charge distribution nearby is

.. math::

   \phi(r) = \frac{1}{4\pi\epsilon_0}
   \int ds\frac{\rho(s)}{|r-s|}.

Tinker uses a variable **electric** (in **chgpot** module) to represent the
the factor :math:`1/(4\pi\epsilon_0)`.
Its default magnitude is 332.063713, which is a constant defined by
**coulomb** (in **units** module), and its units are kcal/mol |ang|/|e2|.
The default value is editable by the *ELECTIRIC* keyword.

.. note::

   Should the value of **coulomb** documented here be out-dated and become
   inconsistent with our code, a pull request will be appreciated.

Expanding :math:`1/|r-s|` in Taylor series, :math:`4\pi\epsilon_0\phi(r)`
can be rewritten as

.. math::
   :label: pot1

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
The third term in :eq:`pot1` can be rewritten as

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
