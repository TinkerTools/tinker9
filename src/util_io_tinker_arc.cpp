#include "mod_box.h"
#include "mod_md.h"
#include "util_array.h"
#include "util_io.h"
#include <ext/tinker/tinker_mod.h>
#include <ext/tinker/tinker_rt.h>
#include <fstream>
#include <sstream>

TINKER_NAMESPACE_BEGIN
void copyin_tinker_arc(const std::string& arcfile, int first1, int last1,
                       int step) {

  if (!(first1 >= 1 && last1 >= first1 && step > 0)) {
    std::string msg = format("Invalid First/Last/Step Values  : {:d}/{:d}/{:d}",
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
      double baxis = std::stod(vs.at(1)); // will throw an exception here if
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
    std::vector<Box> bbuf2;
    if (has_boxsize) {
      bbuf2.resize(tn);
      for (int i = 0; i < tn; ++i) {
        int c = i * 6;
        boxes::xbox = bbuf[c];
        boxes::ybox = bbuf[c + 1];
        boxes::zbox = bbuf[c + 2];
        boxes::alpha = bbuf[c + 3];
        boxes::beta = bbuf[c + 4];
        boxes::gamma = bbuf[c + 5];
        TINKER_RT(lattice)();

        for (int j = 0; j < 3; ++j)
          for (int k = 0; k < 3; ++k) {
            bbuf2[i].lvec[j][k] = boxes::lvec[j][k];
            bbuf2[i].recip[j][k] = boxes::recip[j][k];
          }

        bbuf2[i].volbox = boxes::volbox;
        Box::Shape shape = Box::null;
        if (boxes::orthogonal)
          shape = Box::ortho;
        else if (boxes::monoclinic)
          shape = Box::mono;
        else if (boxes::triclinic)
          shape = Box::tri;
        else if (boxes::octahedron)
          shape = Box::oct;
        bbuf2[i].shape = shape;
      }
      copyin_bytes(trajbox, bbuf2.data(), sizeof(Box) * tn);
    }
    copyin_bytes(trajx, xbuf.data(), sizeof(real) * n * tn);
    copyin_bytes(trajy, ybuf.data(), sizeof(real) * n * tn);
    copyin_bytes(trajz, zbuf.data(), sizeof(real) * n * tn);
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
TINKER_NAMESPACE_END
