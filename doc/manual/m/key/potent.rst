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
An example is as follows:

- BOND |nbsp| A |nbsp| B |nbsp| Force (kcal/mol/|ang2|) |nbsp| Ideal (|ang|)

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

**BOND3 [2 integers & 2 reals]** |not9|

.. index:: BOND3

**BOND4 [2 integers & 2 reals]** |not9|

.. index:: BOND4

**BOND5 [2 integers & 2 reals]** |not9|

.. index:: BOND5

.. seealso::

   :ref:`label-bond`

Angle Bending
-------------

**ANGLETERM [NONE/ONLY]**

.. index:: ANGLETERM

This keyword controls use of the bond angle bending potential energy term.
In the absence of a modifying option, this keyword turns on use of the potential.
The *NONE* option turns off use of this potential energy term.
The *ONLY* option turns off all potential energy terms except for this one.

**ANGLEUNIT [real]**

.. index:: ANGLEUNIT

Sets the scale factor needed to convert the energy value computed by the bond
angle bending potential into units of kcal/mol.
The correct value is force field dependent and typically provided in the header
of the master force field parameter file.
The default value of :math:`(\pi/180)^2` is used, if the *ANGLEUNIT* keyword
is not given in the force field parameter file or the keyfile.

**ANGLE [3 integers & 4 reals]**

.. index:: ANGLE

This keyword provides the values for a single bond angle bending parameter.
The integer modifiers give the atom class numbers for the three kinds of atoms
involved in the angle which is to be defined.
The real number modifiers give the force constant value for the angle and up
to three ideal bond angles in degrees.
In most cases only one ideal bond angle is given, and that value is used for
all occurrences of the specified bond angle.
If all three ideal angles are given, the values apply when the central atom of
the angle is attached to 0, 1 or 2 additional hydrogen atoms, respectively.
This “hydrogen environment” option is provided to implement the
corresponding feature of Allinger’s MM force fields.
The default units for the force constant are kcal/mol/|rad2|, but this can be
controlled via the *ANGLEUNIT* keyword.
An example is as follows:

- ANGLE |nbsp| A1 |nbsp| C |nbsp| A2 |nbsp| Force (kcal/mol/|rad2|) |nbsp| Ideal (deg)

**ANGLEF [3 integers & 3 reals]**

.. index:: ANGLEF

This keyword provides the values for a single bond angle bending parameter
for a SHAPES-style Fourier potential function.
The integer modifiers give the atom class numbers for the three kinds of atoms
involved in the angle which is to be defined.
The real number modifiers give the force constant value for the angle,
the angle shift in degrees, and the periodicity value.
Note that the force constant should be given as the “harmonic” value and not
the native Fourier value.
The default units for the force constant are kcal/mol/|rad2|,
but this can be controlled via the *ANGLEUNIT* keyword.
An example is as follows:

- ANGLEF |nbsp| A1 |nbsp| C |nbsp| A2 |nbsp| Force (kcal/mol/|rad2|) |nbsp| Ideal (deg) |nbsp| Periodicity

**ANGLEP**

.. index:: ANGLEP

**ANGLE-CUBIC [real]**

.. index:: ANGLE-CUBIC

Sets the value (in 1/deg) of the cubic term in the Taylor series expansion
form of the bond angle bending potential energy.
The real number modifier gives the value of the coefficient as a multiple of
the quadratic coefficient.
This term multiplied by the angle bending energy unit conversion factor, the
force constant, and the cube of the deviation of the bond angle from its ideal
value gives the cubic contribution to the angle bending energy.
The default value in the absence of the *ANGLE-CUBIC* keyword is zero;
i.e., the cubic angle bending term is omitted.

**ANGLE-QUARTIC [real]**

.. index:: ANGLE-QUARTIC

Sets the value (in 1/|deg2|) of the quartic term in the Taylor series expansion
form of the bond angle bending potential energy.
The real number modifier gives the value of the coefficient as a multiple of
the quadratic coefficient.
This term multiplied by the angle bending energy unit conversion factor, the
force constant, and the forth power of the deviation of the bond angle from its
ideal value gives the quartic contribution to the angle bending energy.
The default value in the absence of the *ANGLE-QUARTIC* keyword is zero;
i.e., the quartic angle bending term is omitted.

**ANGLE-PENTIC [real]**

.. index:: ANGLE-PENTIC

Sets the value (in 1/|deg3|) of the fifth power term in the Taylor series
expansion form of the bond angle bending potential energy.
The real number modifier gives the value of the coefficient as a multiple of
the quadratic coefficient.
This term multiplied by the angle bending energy unit conversion factor, the
force constant, and the fifth power of the deviation of the bond angle from
its ideal value gives the pentic contribution to the angle bending energy.
The default value in the absence of the *ANGLE-PENTIC* keyword is zero;
i.e., the pentic angle bending term is omitted.

**ANGLE-SEXTIC [real]**

.. index:: ANGLE-SEXTIC

Sets the value (in 1/|deg4|) of the sixth power term in the Taylor series
expansion form of the bond angle bending potential energy.
The real number modifier gives the value of the coefficient as a multiple of
the quadratic coefficient.
This term multiplied by the angle bending energy unit conversion factor, the
force constant, and the sixth power of the deviation of the bond angle from
its ideal value gives the sextic contribution to the angle bending energy.
The default value in the absence of the *ANGLE-SEXTIC* keyword is zero;
i.e., the sextic angle bending term is omitted.

**ANGLE3** |not9|

.. index:: ANGLE3

**ANGLE4** |not9|

.. index:: ANGLE4

**ANGLE5** |not9|

.. index:: ANGLE5

.. seealso::

   :ref:`label-angle`

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
angle parameter as follows:

- IMPROPER |nbsp| D |nbsp| A |nbsp| B |nbsp| C |nbsp| Force (kcal/mol/|rad2|) |nbsp| Ideal (deg)

.. seealso::

   :ref:`label-improp`
