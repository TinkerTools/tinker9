Potential Function Keywords
===========================

.. include:: ../replace.rst

Bond Stretching
---------------

**BONDTYPE [HARMONIC/MORSE]**

.. index:: BONDTYPE

Chooses the functional form of the bond stretching potential.
The *HARMONIC* option selects a Taylor series expansion containing terms from
harmonic through quartic.
The *MORSE* option selects a Morse potential fit to the ideal bond length and
stretching force constant parameter values.
The default is to use the *HARMONIC* potential.

**BONDTERM [NONE/ONLY]**

.. index:: BONDTERM

This keyword controls use of the bond stretching potential energy term.
In the absence of a modifying option, this keyword turns on use of the potential.
The *NONE* option turns off use of this potential energy term.
The *ONLY* option turns off all potential energy terms except for this one.

**BONDUNIT [real]**

.. index:: BONDUNIT

Sets the scale factor needed to convert the energy value computed by the bond
stretching potential into units of kcal/mol. The correct value is force field
dependent and typically provided in the header of the master force field
parameter file. The default value of 1.0 is used, if the *BONDUNIT* keyword
is not given in the force field parameter file or the keyfile.

**BOND [2 integers & 2 reals]**

.. index:: BOND

This keyword provides the values for a single bond stretching parameter.
The integer modifiers give the atom class numbers for the two kinds of atoms
involved in the bond which is to be defined. The real number modifiers give the
force constant value for the bond and the ideal bond length in Angstroms.
An example is as follows.

- BOND |nbsp| A |nbsp| B |nbsp| Force (1/|ang2|) |nbsp| Ideal (|ang|)

**BOND-CUBIC [real]**

.. index:: BOND-CUBIC

Sets the value (in 1/|ang|) of the cubic term in the Taylor series expansion
form of the bond stretching potential energy.
The real number modifier gives the value of
the coefficient as a multiple of the quadratic coefficient.
This term multiplied by the bond stretching energy unit conversion factor,
the force constant, and the cube of the deviation of the bond length from its
ideal value gives the cubic contribution to the bond stretching energy.
The default value in the absence of the *BOND-CUBIC* keyword is zero;
i.e., the cubic bond stretching term is omitted.

**BOND-QUARTIC [real]**

.. index:: BOND-QUARTIC

Sets the value (in 1/|ang2|) of the quartic term in the Taylor series expansion
form of the bond stretching potential energy.
The real number modifier gives the value of
the coefficient as a multiple of the quadratic coefficient.
This term multiplied by the bond stretching energy unit conversion factor,
the force constant, and the forth power of the deviation of the bond length from
its ideal value gives the quartic contribution to the bond stretching energy.
The default value in the absence of the *BOND-QUARTIC* keyword is zero;
i.e., the quartic bond stretching term is omitted.

.. seealso::

   :ref:`label-bond`

Improper Dihedral
-----------------

**IMPROPTERM [NONE/ONLY]**

.. index:: IMPROPTERM

This keyword controls use of the CHARMM-style improper dihedral angle potential
energy term.
In the absence of a modifying option, this keyword turns on use of the potential.
The *NONE* option turns off use of this potential energy term.
The *ONLY* option turns off all potential energy terms except for this one.

**IMPROPUNIT [real]**

.. index:: IMPROPUNIT

Sets the scale factor needed to convert the energy value computed by the
CHARMM-style improper dihedral angle potential into units of kcal/mol.
The correct value is force field dependent and typically provided in the header
of the master force field parameter file. The default value of 1.0 is used,
if the *IMPROPUNIT* keyword is not given in the force field parameter file
or the keyfile.

**IMPROPER [4 integers & 2 reals]**

.. index:: IMPROPER

This keyword provides the values for a single CHARMM-style improper dihedral
angle parameter as follows.

- IMPROPER |nbsp| D |nbsp| A |nbsp| B |nbsp| C |nbsp| Force (1/|rad2|) |nbsp| Ideal (degree)

.. seealso::

   :ref:`label-improp`
