Induced Dipole
==============

Energy
------

:math:`\mu` is the induced dipole in the external field *E*.
The induced field due to the induced dipole is :math:`E^u=-T\mu`, and the
induced dipole is proportional to the total field :math:`E^t`:

.. math::

   \mu = \alpha E^t = \alpha(E+E^u),

where :math:`\alpha` is the polarizability.
Defining :math:`\tilde{T}=\alpha^{-1}+T`, the induced dipole is the solution
to the linear equation

.. math::
   :label: u-dipole

   \tilde{T}\mu = E.

The polarization energy is given by

.. math::
   :label: polar1

   U = -\mu E - \frac{1}{2} \mu E^u + \int_0^\mu d\mu\ \mu E^t_\mu.

On the right-hand side of :eq:`polar1`:

   - the 1st term is the contribution from the external field;
   - the 2nd term is the mutual polarization energy per induced dipole;
   - the 3rd term is the self induced dipole energy and is equal to :math:`\mu\alpha^{-1}\mu/2`.

Finally, the polarization energy is

.. math::
   :label: polar2

   U = -\frac{1}{2}\mu E.

Gradient
--------

With limited numerical precision, the solution to linear equation :eq:`u-dipole`
cannot be fully precise:

.. math::
   :label: u-dipole2

   \tilde{T}\mu = \epsilon + E.

The gradient of the induced dipole can be written in

.. math::

   \mu' = \tilde{T}^{-1}(\epsilon' + E' - \tilde{T}'\mu),

and the polarization gradient is

.. math::

   U' &= -\frac{1}{2} (E\mu' + \mu E') \\
      &= -\frac{1}{2} [\mu\tilde{T}\tilde{T}^{-1}(\epsilon' +E' -\tilde{T}'\mu) +\mu E'] \\
      &= -\frac{1}{2} (\mu\epsilon' +\mu E' -\mu\tilde{T}'\mu +\mu E').

If the convergence of :eq:`u-dipole2` is tight, :math:`\epsilon` and :math:`\epsilon'`
terms will drop. The polarization gradient becomes

.. math::

   U' = -\frac{1}{2} (\mu E' -\mu\tilde{T}'\mu +\mu E').
