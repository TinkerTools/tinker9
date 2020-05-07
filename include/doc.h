/**
 * \defgroup ff  Force Fields and Molecular Mechanics
 *
 *    \defgroup geom  Geometrical Restraints
 *    \ingroup ff
 *
 *    \defgroup vdw  Van der Waals (VDW)
 *    \ingroup ff
 *
 *    \defgroup pme  Particle Mesh Ewald
 *    \ingroup ff
 *    \defgroup charge  Partial Charge Electrostatics
 *    \ingroup ff
 *    \defgroup mpole  Multipole Electrostatics
 *    \ingroup ff
 *    \defgroup polar  AMOEBA Polarization Electrostatics
 *    \ingroup ff
 */


/**
 * \defgroup md  Molecular Dynamics Simulation
 *
 *    \defgroup box  Periodic Boundary Box
 *    \ingroup md
 *    \defgroup spatial  Spatial Decomposition
 *    \ingroup md
 */


/**
 * \defgroup macro  Macros
 */


/**
 * \defgroup test  Unit Tests
 */


/**
 * \defgroup util  Utilities
 *
 *    \defgroup math  Math
 *    \ingroup util
 *    \defgroup rand  Random Number
 *    \ingroup util
 *
 *    \defgroup io  I/O and Text
 *    \ingroup util
 *
 *    \defgroup error  Errors and Exceptions
 *    \ingroup util
 *
 *    \defgroup atomic  Atomic Operations
 *    \ingroup util
 *    \defgroup mem  Memory and Pointers
 *    \ingroup util
 *    \defgroup nvidia  Nvidia GPU
 *    \ingroup util
 *
 *    \defgroup bindc    Fortran Tinker to C Interface
 *    \ingroup util
 */


/**
 * \page keywords  Use of the Keyword Control File
 *
 * Keywords are read from the keyword control file.
 *
 * Several keywords take a list of integer values (atom numbers, for example) as
 * modifiers. For these keywords the integers can simply be listed explicitly
 * and separated by spaces, commas or tabs. If a range of numbers is desired, it
 * can be specified by listing the negative of the first number of the range,
 * followed by a separator and the last number of the range. For example, the
 * keyword line `LIGAND 4 -9 17 23` could be used to add atoms 4, 9 through 17,
 * and 23 to the set of ligand atoms during a Tinker calculation.
 *
 * ## Keywords Grouped by Functionality
 *
 * Listed below are the available Tinker keywords sorted into groups by general
 * function, along with brief descriptions of their actions, possible keyword
 * modifiers, and usage examples.
 *
 * - \subpage kmath  "Mathematical Algorithm Keywords"
 * - \subpage kosrw  "OSRW Keywords"
 */
