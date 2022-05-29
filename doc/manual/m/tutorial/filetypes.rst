Common File Types
=================

You will find *.xyz* and *.key* files under *example* directory and *.prm* files
under *params* directory.

sample.xyz
   The *.xyz* file is the basic Tinker Cartesian coordinates file type.
   It contains a title line followed by one line for each atom in the structure.
   Each line contains: the sequential number within the structure, an atomic
   symbol or name, X-, Y-, and Z-coordinates, the force field atom type number
   of the atom, and a list of the atoms connected to the current atom.
   Except for programs whose basic operation is in torsional space,
   all Tinker calculations are done from some version of the *.xyz* format.

sample.key
   The keyword parameter file always has the extension *.key* and is optionally
   present during Tinker calculations. It contains values for any of a wide
   variety of switches and parameters that are used to change the course of the
   computation from the default. The detailed contents of this file is explained
   in a latter section. If a molecular system specific keyfile, in this case
   *sample.key*, is not present, the the Tinker program will look in the same
   directory for a generic file named *tinker.key*.

sample.prm
   The potential energy parameter files distributed with the Tinker package all
   end in the extension *.prm*, although this is not required by the programs
   themselves. Each of these files contains a definition of the potential energy
   functional forms for that force field as well as values for individual energy
   parameters. For example, the *mm3pro.prm* file contains the energy parameters
   and definitions needed for a protein-specific version of the MM3 force field.
