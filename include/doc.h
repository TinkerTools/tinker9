/**
 * \defgroup bindc    Fortran Tinker to C Interface
 * \defgroup ff       Force Field Terms
 * \defgroup k        Keywords
 * \defgroup macro    Macros
 * \defgroup math     Math
 * \defgroup md       MD Configurations
 * \defgroup nvidia   Nvidia GPU
 * \defgroup spatial  Spatial Decomposition
 * \defgroup test     Unit Tests
 * \defgroup util     Utilities
 */


/**
 * \defgroup elec  Electrostatics
 * \ingroup ff
 *    \defgroup charge  Partial Charge Electrostatics
 *    \ingroup elec
 *    \defgroup mpole   Multipole Electrostatics
 *    \ingroup elec
 *    \defgroup polar   AMOEBA Polarization Electrostatics
 *    \ingroup elec
 *    \defgroup pme     Particle Mesh Ewald
 *    \ingroup elec
 *
 * \defgroup geom  Geometrical Restraints
 * \ingroup ff
 *
 * \defgroup vdw   Van der Waals (VDW)
 * \ingroup ff
 */


/**
 * \defgroup kmath  Mathematical Algorithm Keywords
 * \ingroup k
 * \defgroup kosrw  OSRW Keywords
 * \ingroup k
 */


/**
 * \defgroup atomic Atomic Operations
 * \ingroup util
 * \defgroup box    Periodic Boundary Box
 * \ingroup util
 * \defgroup error  Errors and Exceptions
 * \ingroup util
 * \defgroup io     I/O and Text
 * \ingroup util
 * \defgroup mem    Memory and Pointers
 * \ingroup util
 */


namespace tinker {
/**
 * \brief Integer flags for different kinds of atom data and energy routines.
 */
namespace calc {}


/**
 * \brief Math functions running in parallel.
 */
namespace parallel {}
}
