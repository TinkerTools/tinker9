Potential Function Keywords
===========================

Improper Dihedral Angle
-----------------------

**IMPROPTERM [NONE/ONLY]**

.. index:: IMPROPTERM

This keyword controls use of the CHARMM-style improper dihedral angle potential
energy term. In the absence of a modifying option, this keyword turns on use of
the potential. The *NONE* option turns off use of this potential energy term.
The *ONLY* option turns off all potential energy terms except for this one.

**IMPROPUNIT [real]**

.. index:: IMPROPUNIT

Sets the scale factor needed to convert the energy value computed by the
CHARMM-style improper dihedral angle potential into units of kcal/mole.
The correct value is force field dependent and typically provided in the header
of the master force field parameter file. The default value of 1.0 is used,
if the *IMPROPUNIT* keyword is not given in the force field parameter file
or the keyfile.

**IMPROPER [4 integers & 2 reals]**

.. index:: IMPROPER

This keyword provides the values for a single CHARMM-style improper dihedral
angle parameter.

.. seealso::

   :ref:`label-improp`
