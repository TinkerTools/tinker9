#include "tinkerrt.h"
#include "version.h"
#include <tinker/modcpp.h>
#ifdef _OPENMP
#   include <omp.h>
#endif

namespace tinker {
void initial()
{
   static bool first = true;

   using namespace iounit;
   input = 5;
   iout = 6;

   if (first)
      promo();

   if (first)
      tinker_f_command();
   if (first)
      first = false;

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

   tinker_f_initatom();

   tinker_f_initres();

   using namespace keys;
   nkey = 0;

   using namespace params;
   nprm = 0;

   using namespace atoms;
   n = 0;

   using namespace molcul;
   nmol = 0;

   using namespace cell;
   ncell = 1;

   using namespace align;
   nfit = 0;

   using namespace mutant;
   nmut = 0;

   using namespace zclose;
   nadd = 0;
   ndel = 0;

   using namespace pdb;
   npdb = 0;

   using namespace sequen;
   nseq = 0;
   nchain = 0;

   using namespace files;
   nprior = 0;

   using namespace fft;
   planf = 0;
   planb = 0;

   using namespace inform;
   silent = false;
   verbose = false;
   debug = false;
   abort = false;

   using namespace group;
   use_group = false;

   using namespace bound;
   use_bounds = false;
   use_replica = false;
   use_polymer = false;

   using namespace neigh;
   dovlst = true;
   dodlst = true;
   doclst = true;
   domlst = true;
   doulst = true;

   using namespace bath;
   isothermal = false;
   isobaric = false;

   using namespace virial;
   use_virial = true;

   using namespace rigid;
   use_rigid = false;

   using namespace scales;
   set_scale = false;

   using namespace socket;
   sktstart = false;
   use_socket = false;

   using namespace warp;
   use_smooth = false;
   use_dem = false;
   use_gda = false;
   use_tophat = false;
   use_stophat = false;

   FstrView coordtype = output::coordtype;
   coordtype = "NONE";

   using namespace boxes;
   xbox = 0;
   ybox = 0;
   zbox = 0;
   alpha = 0;
   beta = 0;
   gamma = 0;

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
