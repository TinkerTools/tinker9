#include "gpu/decl_random.h"
#include "gpu/rc.h"
#include "util/fort_str.h"
#include "util/text.h"
#include <chrono>
#include <random>

TINKER_NAMESPACE_BEGIN
namespace gpu {
/// @return  Zero if Tinker key file does not consist of a RANDOMSEED keyword
static int read_tinker_randomseed_() {
  int seed = 0;
  for (int i = 0; i < keys::nkey; ++i) {
    fstr_view record = keys::keyline[i];
    auto vs = Text::split(record.trim());
    if (vs.size()) {
      std::string keyword = vs.at(0);
      Text::upcase(keyword);
      if (keyword == "RANDOMSEED" && vs.size() > 1) {
        seed = std::stoi(vs.at(1));
        seed = std::max(1, seed);
      }
    }
  }
  return seed;
}

static std::default_random_engine generator_;
void random_data(rc_t rc) {
  if (rc & rc_copyin) {
    int seed = read_tinker_randomseed_();
    if (seed == 0) {
      auto now = std::chrono::system_clock::now();
      auto tt = std::chrono::system_clock::to_time_t(now);
      auto local_tm = std::localtime(&tt);
      int year = local_tm->tm_year % 10;
      int month = local_tm->tm_mon;
      int day = local_tm->tm_mday;
      int hour = local_tm->tm_hour;
      int minute = local_tm->tm_min;
      int second = local_tm->tm_sec;
      seed = 32140800 * year + 2678400 * (month - 1);
      seed += 86400 * (day - 1) + 3600 * hour;
      seed += 60 * minute + second;
    }
    generator_.seed(seed);
  }
}

static std::uniform_real_distribution<double> uniformd_(0, 1);
static std::uniform_real_distribution<float> uniformf_(0, 1);
double random_double() { return uniformd_(generator_); }
float random_float() { return uniformf_(generator_); }

static std::normal_distribution<double> normald_(0, 1);
static std::normal_distribution<float> normalf_(0, 1);
double normal_double() { return normald_(generator_); }
float normal_float() { return normalf_(generator_); }

double chi_squared_double(int k) {
  std::chi_squared_distribution<double> chi_squared_(k);
  return chi_squared_(generator_);
}
float chi_squared_float(int k) {
  std::chi_squared_distribution<float> chi_squared_(k);
  return chi_squared_(generator_);
}
}
TINKER_NAMESPACE_END
