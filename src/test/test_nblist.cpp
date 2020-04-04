#include "error.h"
#include "files.h"
#include "io_print.h"
#include "io_text.h"
#include "md.h"
#include "nblist.h"
#include "platform.h"
#include "rc_man.h"
#include "test.h"
#include "test_rt.h"
#include "tinker_rt.h"
#include <fstream>
#include <set>
#include <sstream>
#include <tinker/detail/boxes.hh>


using namespace TINKER_NAMESPACE;


// static const int ans[5][216][70];
#include "test_nblist.hh"


namespace {
// if all integers in ref can be found in array
bool find_match(const int* array, int na, int iframe, int iatom)
{
   std::set<int> array_set(array, array + na);
   std::vector<int> missing;

   for (int i = 0; i < nmax; ++i) {
      int r = ans[iframe][iatom][i];
      if (r >= 0) {
         if (!array_set.count(r)) {
            missing.push_back(r);
         }
      }
   }

   if (missing.size() == 0) {
      return true;
   } else {
      std::string line =
         format("The following neighbors of atom %3d in frame %d are missing:",
                iatom + 1, iframe + 1);
      for (auto i : missing) {
         line += format(" %3d", i + 1);
      }
      line += "\n";
      print(stdout, line);
      return false;
   }
}


void goto_frame(int idx0)
{
   assert(calc::traj & rc_flag);
   x = trajx + n * idx0;
   y = trajy + n * idx0;
   z = trajz + n * idx0;
   const Box& b = trajbox[idx0];
   set_default_box(b);
}


void copyin_arc_file(const std::string& arcfile, int first1, int last1,
                     int step)
{

   if (!(first1 >= 1 && last1 >= first1 && step > 0)) {
      std::string msg = format("Invalid First/Last/Step Values  : %d/%d/%d",
                               first1, last1, step);
      TINKER_THROW(msg);
   }

   // (1, 7, 2) -> 4
   // (1, 8, 2) -> 4
   int tn = (last1 - first1) / step + 1;
   assert(tn <= trajn);

   double xyzsave[3], anglesave[3];
   xyzsave[0] = boxes::xbox;
   xyzsave[1] = boxes::ybox;
   xyzsave[2] = boxes::zbox;
   anglesave[0] = boxes::alpha;
   anglesave[1] = boxes::beta;
   anglesave[2] = boxes::gamma;

   std::ifstream iarc(arcfile);
   if (iarc) {
      std::string line;
      std::getline(iarc, line); // n and title

      std::getline(iarc, line); // either box size or first atom
      bool has_boxsize;
      try {
         auto vs = Text::split(line);
         std::stod(vs.at(1)); // will throw an exception here if
                              // this file does not include box size
         has_boxsize = true;
      } catch (std::invalid_argument&) {
         has_boxsize = false;
      }
      int nlines_to_skip = n + 1;
      if (has_boxsize)
         nlines_to_skip += 1;

      std::vector<double> bbuf;
      if (has_boxsize)
         bbuf.resize(tn * 6);
      std::vector<real> xbuf(tn * n), ybuf(tn * n), zbuf(tn * n);

      auto skip = [](std::ifstream& f, int nl) {
         for (int i = 0; i < nl; ++i)
            f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      };

      // rewind
      iarc.clear();
      iarc.seekg(0);

      int current = 0;
      int dummy1;
      char dummy2[8];
      for (int i = 1; i <= last1; i += step) {
         bool require1 = (i >= first1);
         bool require2 = ((i - first1) % step == 0);
         if (require1 && require2) {
            skip(iarc, 1); // skip n and title
            if (has_boxsize) {
               std::getline(iarc, line);
               std::istringstream ss(line);
               int c = 6 * current;
               ss >> bbuf[c] >> bbuf[c + 1] >> bbuf[c + 2] >> bbuf[c + 3] >>
                  bbuf[c + 4] >> bbuf[c + 5];
            }
            int off = current * n;
            for (int j = 0; j < n; ++j) {
               std::getline(iarc, line);
               std::istringstream ss(line);
               ss >> dummy1 >> dummy2 >> xbuf[off + j] >> ybuf[off + j] >>
                  zbuf[off + j];
            }
            ++current;
         } else {
            skip(iarc, nlines_to_skip);
         }
      }

      // copyin
      if (has_boxsize) {
         for (int i = 0; i < tn; ++i) {
            int c = i * 6;
            boxes::xbox = bbuf[c];
            boxes::ybox = bbuf[c + 1];
            boxes::zbox = bbuf[c + 2];
            boxes::alpha = bbuf[c + 3];
            boxes::beta = bbuf[c + 4];
            boxes::gamma = bbuf[c + 5];
            TINKER_RT(lattice)();
            get_tinker_box_module(trajbox[i]);
         }
      }
      darray::copyin(PROCEED_NEW_Q, n * tn, trajx, xbuf.data());
      darray::copyin(PROCEED_NEW_Q, n * tn, trajy, ybuf.data());
      darray::copyin(WAIT_NEW_Q, n * tn, trajz, zbuf.data());
   } else {
      std::string msg = "Cannot Open File ";
      msg += arcfile;
      TINKER_THROW(msg);
   }

   boxes::xbox = xyzsave[0];
   boxes::ybox = xyzsave[1];
   boxes::zbox = xyzsave[2];
   boxes::alpha = anglesave[0];
   boxes::beta = anglesave[1];
   boxes::gamma = anglesave[2];
}
}


TEST_CASE("NBList-ArBox", "[ff][nblist][arbox]")
{
   int usage_ = calc::xyz | calc::mass | calc::traj | calc::energy;

   const char* k = "test_arbox.key";
   const char* x1 = "test_arbox.arc";
   const char* p = "amoeba09.prm";

   std::string k0 = arbox_key;
   TestFile fke(k, k0);

   TestFile fx1(x1, arbox_arc);
   TestFile fpr(p, commit_6fe8e913::amoeba09_prm);

   const char* argv[] = {"dummy", x1};
   int argc = 2;
   test_begin_with_args(argc, argv);

   rc_flag = usage_;
   pltfm_config = ACC_PLTFM; // to always use nblist
   trajn = 5;
   initialize();

   copyin_arc_file(x1, 1, 5, 1);

   const int maxnlst = vlist_unit->maxnlst;
   std::vector<int> nlst;
   nlst.resize(n);
   std::vector<int> lst;
   lst.resize(n * maxnlst);

   for (int ifr = 0;;) {
      darray::copyout(PROCEED_NEW_Q, n, nlst.data(), vlist_unit->nlst);
      darray::copyout(WAIT_NEW_Q, n * maxnlst, lst.data(), vlist_unit->lst);

      for (int iatom = 0; iatom < n; ++iatom)
         REQUIRE(find_match(&lst[iatom * maxnlst], nlst[iatom], ifr, iatom));

      ++ifr;
      // update nblist
      if (ifr < trajn) {
         goto_frame(ifr);
         refresh_neighbors();
      } else
         break;
   }

   finish();
   test_end();
}
