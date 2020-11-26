Potential Energy Functions
==========================

.. _label-bond:

Bond Stretching
---------------

Bond term is an empirical function of the deviation of the bond length from
its ideal value (:math:`\Delta b = b - b_0`).
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

Parameter *a* is hardwired to 2. Following equation
:math:`a = \sqrt{\frac{k}{2 D_e}}` and the Tinker convention to include 1/2 in
the force constant, *De* is hardwired to *k*/4.

.. _label-angle:

Angle Bending
-------------

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

It is worth noting that there is no 1/2 coefficient before the force
coefficient, which may be different in other software packages.
