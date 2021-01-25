#include "promo.h"
#include "subroutine.h"
#include "tool/io_fort_str.h"
#include <tinker/tinker_mod.h>
#ifdef _OPENMP
#   include <omp.h>
#endif


extern "C"
{
   void command_();
   void initatom_();
   void initres_();
}


namespace tinker {
void initial()
{
   static bool first = true;


   // default unit numbers for input and output
   using namespace iounit;
   input = 5;
   iout = 6;


   // display program banner and copyright notice
   if (first)
      promo();


   // command line arguments to the program
   if (first)
      command_();
   if (first)
      first = false;


   // cores, thread count and options for OpenMP
   using namespace openmp;
   nproc = 1;
   nthread = 1;
#ifdef _OPENMP
   nproc = omp_get_num_procs();
   nthread = nproc;
   omp_set_num_threads(nthread);
   omp_set_nested(true);
#   ifdef TINKER_ICPC
   // Intel compiler extensions to OpenMP standard,
   // 268435456 bytes is 2**28 bytes, or 256 MB
   kmp_set_stacksize_s(268435456);
   kmp_set_blocktime(0);
#   endif
#endif


   // atomic symbols for elements
   initatom_();


   // names of biopolymer residue types
   initres_();


   // number of lines in the keyfile
   using namespace keys;
   nkey = 0;


   // number of lines in the parameter file
   using namespace params;
   nprm = 0;


   // number of atoms in the system
   using namespace atoms;
   n = 0;


   // number of molecules in the system
   using namespace molcul;
   nmol = 0;


   // number of unit cell and replicates
   using namespace cell;
   ncell = 1;


   // number of atoms used in superposition
   using namespace align;
   nfit = 0;


   // number of mutated atoms in the system
   using namespace mutant;
   nmut = 0;


   // number of bonds added or deleted from Z-matrix
   using namespace zclose;
   nadd = 0;
   ndel = 0;


   // number of atoms in Protein Data Bank format
   using namespace pdb;
   npdb = 0;


   // number of residues and chains in biopolymer sequence
   using namespace sequen;
   nseq = 0;
   nchain = 0;


   // highest numbered previous cycle file
   using namespace files;
   nprior = 0;


   // pointer initialization for FFTW plans
   using namespace fft;
   planf = 0;
   planb = 0;


   // flags for information levels within the program
   using namespace inform;
   silent = false;
   verbose = false;
   debug = false;
   abort = false;


   // flag for use of atom groups
   using namespace group;
   use_group = false;


   // flags for use of periodic boundaries
   using namespace bound;
   use_bounds = false;
   use_replica = false;
   use_polymer = false;


   // flags for rebuilding of neighbor lists
   using namespace neigh;
   dovlst = true;
   dodlst = true;
   doclst = true;
   domlst = true;
   doulst = true;


   // flags for temperature and pressure baths
   using namespace bath;
   isothermal = false;
   isobaric = false;


   // flag for use of internal virial
   using namespace virial;
   use_virial = true;


   // flag for use of rigid bodies
   using namespace rigid;
   use_rigid = false;


   // flag to show setting of optimization scale factors
   using namespace scales;
   set_scale = false;


   // flags for external Java socket communication
   using namespace socket;
   sktstart = false;
   use_socket = false;


   // flags for potential energy smoothing
   using namespace warp;
   use_smooth = false;
   use_dem = false;
   use_gda = false;
   use_tophat = false;
   use_stophat = false;


   // type of coordinates file
   fstr_view coordtype = output::coordtype;
   coordtype = "NONE";


   // default values for unit cell dimensions
   using namespace boxes;
   xbox = 0;
   ybox = 0;
   zbox = 0;
   alpha = 0;
   beta = 0;
   gamma = 0;


   // default values used by optimizations
   using namespace linmin;
   using namespace minima;
   fctmin = 0;
   maxiter = 0;
   nextiter = 0;
   iprint = -1;
   iwrite = -1;
   stpmax = 0;
}
}
