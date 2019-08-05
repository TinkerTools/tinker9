#include "random.h"
#include "rc_man.h"
#include "util_io.h"
#include <ext/tinker/tinker_mod.h>

#include <chrono>
#include <random>

TINKER_NAMESPACE_BEGIN
/// @return
/// Zero if Tinker key file does not contain a RANDOMSEED keyword
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
void random_data(rc_op op) {
  if (op & rc_init) {
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

static std::gamma_distribution<double> gammad_(1, 1);
static std::gamma_distribution<float> gammaf_(1, 1);
double chi_squared_double(int k) {
  gammad_.param(std::gamma_distribution<double>::param_type(0.5 * k, 2));
  return gammad_(generator_);
}
float chi_squared_float(int k) {
  gammaf_.param(std::gamma_distribution<float>::param_type(0.5f * k, 2));
  return gammaf_(generator_);
}
TINKER_NAMESPACE_END
