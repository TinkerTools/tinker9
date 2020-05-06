#include "io_print.h"
#include "subroutine.h"
#include <tinker/detail/inform.hh>


extern "C"
{
   void attach_();
   void active_();

   void bonds_();
   void angles_();
   void torsions_();
   void bitors_();
   void rings_();

   void field_();

   void unitcell_();
   void lattice_();
   void polymer_();
   void cutoffs_();

   void flatten_();

   void katom_();

   void molecule_();
   void cluster_();

   void orbital_();

   void kbond_();
   void kangle_();
   void kstrbnd_();
   void kurey_();
   void kangang_();

   void kopbend_();
   void kopdist_();
   void kimprop_();
   void kimptor_();

   void ktors_();
   void kpitors_();
   void kstrtor_();
   void kangtor_();
   void ktortor_();

   void kvdw_();
   void krepel_();
   void kdisp_();

   void kcharge_();
   void kdipole_();
   void kmpole_();
   void kpolar_();
   void kchgtrn_();

   void ksolv_();
   void kmetal_();
   void korbit_();
   void kgeom_();
   void kextra_();

   void kewald_();

   void shakeup_();

   void mutate_();

   void fatal_();
}


namespace tinker {
void mechanic()
{
   // set the bonded connectivity lists and active atoms
   attach_();
   active_();


   // find bonds, angles, torsions, bitorsions and small rings
   bonds_();
   angles_();
   torsions_();
   bitors_();
   rings_();


   // get the base force field from parameter file and keyfile
   field_();


   // find unit cell type, lattice parameters and cutoff values
   unitcell_();
   lattice_();
   polymer_();
   cutoffs_();


   // setup needed for potential energy smoothing methods
   flatten_();


   // assign atom types, classes and other atomic information
   katom_();


   // assign atoms to molecules and set the atom groups
   molecule_();
   cluster_();


   // find any pisystem atoms, bonds and torsional angles
   orbital_();


   // assign bond, angle and cross term potential parameters
   kbond_();
   kangle_();
   kstrbnd_();
   kurey_();
   kangang_();


   // assign out-of-plane deformation potential parameters
   kopbend_();
   kopdist_();
   kimprop_();
   kimptor_();


   // assign torsion and torsion cross term potential parameters
   ktors_();
   kpitors_();
   kstrtor_();
   kangtor_();
   ktortor_();


   // assign van der Waals, repulsion and dispersion parameters
   kvdw_();
   krepel_();
   kdisp_();


   // assign electrostatic interaction potential parameters
   kcharge_();
   kdipole_();
   kmpole_();
   kpolar_();
   kchgtrn_();


   // assign solvation, metal, pisystem and restraint parameters
   ksolv_();
   kmetal_();
   korbit_();
   kgeom_();
   kextra_();


   // assign electrostatic and dispersion Ewald sum parameters
   kewald_();


   // set any holonomic interatomic distance constraints
   shakeup_();


   // set hybrid parameter values for free energy perturbation
   mutate_();


   // quit if essential parameter information is missing
   if (inform::abort) {
      print(stdout,
            "\n MECHANIC  --  Some Required Potential Energy"
            " Parameters are Undefined\n");
      fatal_();
   }
}
}
