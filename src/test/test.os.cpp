#include "test/test.os.h"
#include <cstdio>
#include <fstream>

TINKER_NAMESPACE_BEGIN
namespace test {
file_gen::file_gen(const std::string& name, const std::string& content)
    : good_(false), name_(name) {
  std::ofstream fout(name);
  good_ = fout.is_open();
  if (good_) {
    fout << content;
  }
}

file_gen::~file_gen() {
  if (good_)
    std::remove(name_.c_str());
}
}
TINKER_NAMESPACE_END
