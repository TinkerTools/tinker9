Potential Energy Functions
==========================

.. include:: ../replace.rst

.. _label-bond:

Bond Stretching
---------------

Bond term is an empirical function of bond deviating from the ideal
bond length (:math:`\Delta b = b - b_0`).
To support the AMOEBA force field model, Tinker includes the 3rd and 4th order
terms.

.. math::

   U = k\Delta b^2(1 + k_3\Delta b + k_4\Delta b^2).

Setting 3rd and 4th order coefficients to zero will give the harmonic
functional form.

.. note::

   Different from Hooke's Law (:math:`U = k x^2/2`), Tinker usually drops
   the coefficient 1/2.

The Morse oscillator is also implemented in Tinker:

.. math::

   U = D_e [1 - \exp(-a\Delta b)]^2.

Parameter *a* is hardwired to 2 by approximation. Following equation
:math:`a = \sqrt{\frac{k}{2 D_e}}` and the Tinker convention to include 1/2 in
the force constant, *De* is *k*/4.

.. _label-angle:

Angle Bending
-------------

Similar to bond stretching, angle bending term is also an empirical
function of angle deviating from the ideal angle value
(:math:`\Delta\theta=\theta-\theta_0`).
Terms from cubic to sextic are added to generalize the *HARMONIC* functional from.

.. math::

   U = k\Delta\theta^2(1+k_3\Delta\theta+k_4\Delta\theta^2
                        +k_5\Delta\theta^3+k_6\Delta\theta^4).

MMFF force field has a special treatment for *LINEAR* angle,
e.g., carbon dioxide.
Since the ideal angle should always be :math:`\pi` rad, the deviation can be
approximated by

.. math::

   \Delta\theta=\theta-\pi=2(\frac{\theta}{2}-\frac{\pi}{2})\sim
   2\sin(\frac{\theta}{2}-\frac{\pi}{2})=-2\cos\frac{\theta}{2}.

Only keeping the quadratic term, the angle bending term can be simplified to

.. math::

   U = 2k(1+\cos\theta).

The *LINEAR* angle type is a special case of the SHAPES-style Fourier potential
function :cite:`shapes-ff` with 1-fold periodicity, which is referred to as the
*FOURIER* angle type in Tinker jargon and has the following form

.. math::

   U = 2k(1+\cos(n\theta-\theta_0)).

.. _label-strbnd:

Stretch-Bend Coupling
---------------------

The functional forms for bond stretching, angle bending, and stretch-bend
coupling are those of the MM3 force field :cite:`mm3-ff`:

.. math::

   U = (k_1\Delta b_1 + k_2\Delta b_2)\Delta\theta.

Even though force constants *k1* and *k2* are implemented as two independent
variables in Tinker, they were treated as the same variable in the original
literature.

.. _label-improp:

Improper Dihedral
-----------------

Commonly used in the CHARMM force fields, this potential function is meant to
keep atoms planar. The ideal angle :math:`\varphi_0` defined by dihedral
*D-A-B-C* will always be zero degrees. *D* is the trigonal atom, *A-B-C* are the
peripheral atoms. In the original CHARMM parameter files, the trigonal atom is
often listed last as *C-B-A-D*.

Some of the improper angles are "double counted" in the CHARMM protein
parameter set. Since Tinker uses only one improper parameter per site, we have
doubled these force constants in the Tinker version of the CHARMM parameters.
Symmetric parameters, which are the origin of the "double counted" CHARMM
values, are handled in the Tinker package by assigning all symmetric states and
using the Tinker force constant divided by the symmetry number.

The harmonic functional form implemented in Tinker is

.. math::

   U = k(\varphi-\varphi_0)^2.

There is no coefficient 1/2 before the force coefficient,
which may be different in other software packages.
