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

Expanding :math:`1/|r-s|` in Taylor series, the potential can be rewritten as
