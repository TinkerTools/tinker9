#include "files.h"
#include "mathfunc.h"
#include "md.h"
#include "test.h"
#include "test_rt.h"
using namespace TINKER_NAMESPACE;


TEST_CASE("CUDAReduce", "[util][math][reduce]")
{
   rc_flag = calc::xyz | calc::vmask;


   const int N = 500;
   const int H = 6;
   const int D = 8;
   std::vector<int> vi(N, 0), vi2(N * D, 0);
   std::vector<float> vf(N, 0), vf2(N * D, 0);
   std::vector<double> vd(N, 0), vd2(N * D, 0);
   std::vector<unsigned long long> vu(N, 0), vu2(N * D, 0);
   for (int i = 0; i < N; ++i) {
      int j = i + 1;
      vi[i] = j;
      vf[i] = j;
      vd[i] = j;
      vu[i] = j;
      for (int k = 0; k < H; ++k) {
         int ii = i * D + k;
         int jj = j * (k + 1);
         vi2[ii] = jj;
         vf2[ii] = jj;
         vd2[ii] = jj;
         vu2[ii] = jj;
      }
   }
   int* di;
   float* df;
   double* dd;
   unsigned long long* du;
   float(*df2)[D];
   double(*dd2)[D];
   unsigned long long(*du2)[D];


   int refi = 125250;
   float reff = refi;
   double refd = refi;
   unsigned long long refu = refi;
   float reff2[H] = {reff * 1, reff * 2, reff * 3,
                     reff * 4, reff * 5, reff * 6};
   double refd2[H] = {refd * 1, refd * 2, refd * 3,
                      refd * 4, refd * 5, refd * 6};
   unsigned long long refu2[H] = {refu * 1, refu * 2, refu * 3,
                                  refu * 4, refu * 5, refu * 6};
   int ai;
   float af;
   double ad;
   unsigned long long au;
   float af2[H];
   double ad2[H];
   unsigned long long au2[H];


   const char* k = "test_trpcage.key";
   const char* x1 = "test_trpcage.xyz";
   const char* p = "amoebapro13.prm";
   std::string k0 = trpcage_key;
   k0 += "\nbondterm only\n";
   k0 += "\ngpu-package cuda\n";
   TestFile fke(k, k0);
   TestFile fx1(x1, trpcage_xyz);
   TestFile fpr(p, commit_6fe8e913::amoebapro13_prm);
   const char* argv[] = {"dummy", x1};
   int argc = 2;
   test_begin_with_args(argc, argv);
   initialize();


   device_array::allocate(N, &di, &df, &dd, &du);
   device_array::allocate(N, &df2, &dd2, &du2);
   device_array::copyin(PROCEED_NEW_Q, N, di, vi.data());
   device_array::copyin(PROCEED_NEW_Q, N, df, vf.data());
   device_array::copyin(PROCEED_NEW_Q, N, dd, vd.data());
   device_array::copyin(PROCEED_NEW_Q, N, du, vu.data());
   device_array::copyin(PROCEED_NEW_Q, N, df2, vf2.data());
   device_array::copyin(PROCEED_NEW_Q, N, dd2, vd2.data());
   device_array::copyin(WAIT_NEW_Q, N, du2, vu2.data());


   ai = parallel::reduce_sum(di, N, false);
   REQUIRE(ai == refi);

   af = parallel::reduce_sum(df, N, false);
   REQUIRE(af == reff);


   ad = parallel::reduce_sum(dd, N, false);
   REQUIRE(ad == refd);


   au = parallel::reduce_sum(du, N, false);
   REQUIRE(au == refu);


   parallel::reduce_sum2(af2, df2, N, false);
   for (int j = 0; j < H; ++j)
      REQUIRE(af2[j] == reff2[j]);


   parallel::reduce_sum2(ad2, dd2, N, false);
   for (int j = 0; j < H; ++j)
      REQUIRE(ad2[j] == refd2[j]);


   parallel::reduce_sum2(au2, du2, N, false);
   for (int j = 0; j < H; ++j)
      REQUIRE(au2[j] == refu2[j]);


   device_array::deallocate(di, df, dd, du);
   device_array::deallocate(df2, dd2, du2);


   finish();
   test_end();
}
