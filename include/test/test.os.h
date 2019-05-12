#ifndef TINKER_TEST_TEST_OS_H_
#define TINKER_TEST_TEST_OS_H_

#include "util/fort.rt.h"
#include <string>

TINKER_NAMESPACE_BEGIN
namespace test {
class file_gen {
private:
  bool good_;
  std::string name_;

public:
  file_gen(const std::string& name, const std::string& content);
  ~file_gen();
};

void test_begin_1_xyz(int argc, const char** argv);
void test_end();
}
TINKER_NAMESPACE_END

#endif
