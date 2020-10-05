Integrators
===========

.. _label-verlet:

Verlet Integrator
-----------------

.. _label-respa:

RESPA Integrator
----------------

.. _label-nose-hoover:

Extended Nosé-Hoover Chain by MTK
---------------------------------

.. index:: NOSE-HOOVER

Several methods for NVT and NPT ensembles were discussed in MTK [#Martyna1996]_.

======  ===============  ======
Number  Sections in MTK  Method
======  ===============  ======
1a      2.1 4.3          NVT
2a      2.2 4.4          NPT (isotropic cell fluctuations)
3a      2.3 4.5          NPT (full cell fluctuations)
4a      5.2              XO-RESPA
4b      5.2              XI-RESPA
1b      5.3              RESPA 1a
2b      5.4              RESPA 2a
3b      5.4              RESPA 3a
======  ===============  ======

The isothermal-isobaric integrator implemented in Fortran Tinker and here is
NPT-XO (#2a-4a).

.. tip::

   MTK Nosé-Hoover Chain can be enabled by keywords

   .. code-block:: text

      integrator nose-hoover

   or

   .. code-block:: text

      thermostat nose-hoover
      barostat   nose-hoover

   with the NPT option in the `dynamic` program.

.. _label-lpiston:

Langevin Piston
---------------

The Langevin piston method for constant pressure [#Feller1995]_ is
integrated in the Leapfrog framework.

.. tip::

   Langevin Piston can be enabled by keywords

      .. code-block:: text

         integrator lpiston

   or

      .. code-block:: text

         thermostat lpiston
         barostat   lpiston

   with the NPT option in the `dynamic` program.

.. [#Martyna1996]
   Martyna, G. J.; Tuckerman, M. E.; Tobias, D. J. and Klein, M. L.
   `Mol. Phys., 87, 1117 (1996) <https://doi.org/10.1080/00268979600100761>`_

.. [#Feller1995]
   Feller, S. E.; Zhang, Y.; Pastor, R. W. and Brooks, B. R.
   `J. Chem. Phys., 103, 4613 (1995) <https://doi.org/10.1063/1.470648>`_
