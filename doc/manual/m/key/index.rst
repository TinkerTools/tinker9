Use of the Keyword Control File
===============================

Keywords are read from the keyword control file.

Several keywords take a list of integer values (atom numbers, for example) as
modifiers. For these keywords the integers can simply be listed explicitly
and separated by spaces, commas or tabs. If a range of numbers is desired, it
can be specified by listing the negative of the first number of the range,
followed by a separator and the last number of the range. For example, the
keyword line **LIGAND 4 -9 17 23** could be used to add atoms 4, 9 through 17,
and 23 to the set of ligand atoms during a Tinker calculation.

Below are the available Tinker keywords sorted into groups by general
function, along with brief descriptions of their actions, possible keyword
modifiers, and usage examples.

.. toctree::

   fep
   math
