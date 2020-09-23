Extended Nosé-Hoover Chain by MTK
=================================

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

.. [#Martyna1996]
   Martyna, G. J.; Tuckerman, M. E.; Tobias, D. J. and Klein, M. L.
   `Mol. Phys., 87, 1117 (1996) <https://doi.org/10.1080/00268979600100761>`_
