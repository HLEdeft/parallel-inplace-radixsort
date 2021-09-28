

// #include "random_shuffle.h"
// #include "seq.h"
// #include "sequence.h"
#include "get_time.h"
#include "sequenceIO.h"
#include "radixSort.h"

#include "parseCommandLine.h"

#include "ska_sort.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <deque>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <unordered_map>

using namespace std;
using namespace benchIO;

typedef std::chrono::high_resolution_clock high_res_clock;
template <typename K, typename V> struct mypair {
  K first;
  V second;
};
typedef pair<uintT, uintT> uintTPair;
int uniform[] = {0,      0,       0,        0,         10,        100,   1000,
                 5000,   7000,    8000,     10000,     15000,     20000, 50000,
                 100000, 1000000, 10000000, 100000000, 1000000000};
int zipfan[] = {10000, 100000, 1000000, 10000000, 100000000, 1000000000};
double exp_lambda[]={1,0.001,0.0003,0.0002,0.00015,0.0001,0.00001};
inline uint64_t hash64_2(uint64_t x) {
  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  x = x ^ (x >> 31);
  return x;
}
template <class F> void check_sorted(uintT *array, long length, F f) {
  uintT prev = array[0];
  parallel_for(long i = 1; i < length; i++) {
    if (array[i] < array[i - 1]) {
      printf("FAIL: ARRAY NOT SORTED!\n");
      exit(0);
    }
  }
  printf("array is sorted!\n");
}

template <class T, class F>
void check_sorted_pairs(pair<uintT, T> *array, long length, F f) {
  bool sorted = true;
  long stability_error = 0;
  parallel_for(long i = 1; i < length; i++) {
    if (f(array[i]) < f(array[i - 1])) {
      sorted = false;
    }
  }
  if (!sorted) {
    printf("FAIL: ARRAY NOT SORTED!\n");
  } else {
    printf("array is sorted!\n");
  }
}
void scan_inplace__(uintT *in, uintT n) {
  if (n <= 10000) {
    // cout << "THRESHOLDS: " << THRESHOLDS << endl;
    for (size_t i = 1; i < n; i++) {
      in[i] += in[i - 1];
    }
    return;
  }
  uintT root_n = (uintT)sqrt(n); // split the array into root n blocks
  uintT *offset = new uintT[root_n - 1];

  // parallel_for (0, root_n-1, [&] (size_t i) {
  parallel_for(size_t i = 0; i < root_n - 1; i++) {
    offset[i] = 0;
    for (size_t j = i * root_n; j < (i + 1) * root_n; j++) {
      offset[i] += in[j];
    }
  }

  for (size_t i = 1; i < root_n - 1; i++)
    offset[i] += offset[i - 1];

  // prefix sum for each subarray
  // parallel_for (0, root_n, [&] (size_t i) {
  parallel_for(size_t i = 0; i < root_n; i++) {
    if (i == root_n - 1) { // the last one
      for (size_t j = i * root_n + 1; j < n; j++) {
        in[j] += in[j - 1];
      }
    } else {
      for (size_t j = i * root_n + 1; j < (i + 1) * root_n; j++) {
        in[j] += in[j - 1];
      }
    }
  }

  // parallel_for (1, root_n, [&] (size_t i) {
  parallel_for(size_t i = 1; i < root_n; i++) {
    if (i == root_n - 1) {
      for (size_t j = i * root_n; j < n; j++) {
        in[j] += offset[i - 1];
      }
    } else {
      for (size_t j = i * root_n; j < (i + 1) * root_n; j++) {
        in[j] += offset[i - 1];
      }
    }
  }

  delete[] offset;
}
void zipfian_generator_int64_( // mypair<uintT, uintT> *A,
    uintT zipf_s, int line) {
  // pbbs::sequence<uintT> nums(zipf_s); // in total zipf_s kinds of keys
  int n = 1000000000;
  std::vector<std::pair<uintT, uintT>> A(n);
  uintT *nums = new uintT[zipf_s];

  /* 1. making nums[] array */
  uintT number = (uintT)(n / log(n)); // number= n/ln(n)
  parallel_for(uintT i = 0; i < zipf_s; i++) {
    nums[i] = (uintT)(number / (i + 1));
  }

  // the last nums[zipf_s-1] should be (n - \sum{zipf_s-1}nums[])
  uintT offset = 0;
  for (uintT i = 0; i < zipf_s; i++) {
    offset += nums[i];
  }
  // uintT offset = pbbs::reduce(nums, pbbs::addm<uintT>()); // cout << "offset
  // = " << offset << endl;
  nums[0] += (n - offset);

  uintT offset_ = 0;
  for (uintT i = 0; i < zipf_s; i++) {
    offset_ += nums[i];
  }
  // checking if the sum of nums[] equals to n
  if ((uintT)offset_ == (uintT)n) {
    cout << "sum of nums[] == n" << endl;
  }

  /* 2. scan to calculate position */
  uintT *addr = new uintT[zipf_s];
  // parallel_for (0, zipf_s, [&] (uintT i) {
  parallel_for(uintT i = 0; i < zipf_s; i++) { addr[i] = nums[i]; }
  scan_inplace__(addr, zipf_s); // store all addresses into addr[]

  /* 3. distribute random numbers into A[i].first */
  for (size_t i = 0; i < zipf_s; i++) {
    size_t st = (i == 0) ? 0 : addr[i - 1],
           ed = (i == zipf_s - 1) ? n : addr[i];
    for (size_t j = st; j < ed; j++) {
      A[j].first = hash64_2(i);
    }
  }
  // parallel_for (0, n, [&] (uintT i){
  parallel_for(size_t i = 0; i < n; i++) { A[i].second = hash64_2(i); }

  /* 4. shuffle the keys */
  random_shuffle(A.begin(), A.end());

  delete[] nums;
  delete[] addr;

  uintT length = n;
  uintTPair *control_array = newA(uintTPair, length);
  parallel_for(uintT i = 0; i < length; i++) {
    control_array[i].first = A[i].first;
    control_array[i].second = A[i].second;
  }
  uintTPair *array = newA(uintTPair, length);
  std::cout << "runing test of zipfan input with line " << line << std::endl;
  for (uintT round = 0; round < 5; round++) {
    parallel_for(uintT i = 0; i < length; i++) {
      array[i].first = control_array[i].first;
      array[i].second = control_array[i].second;
    }
    timer t;
    t.start();
    // auto start = high_res_clock::now();
    parallelIntegerSort(array, length, utils::firstF<uintT, uintT>());
    t.stop();
    cout << t.get_total() << endl;
    check_sorted_pairs(array, length, utils::firstF<uintT, uintT>());
    // check_sorted_pairs(array, length, utils::firstF<uintT, uintT>());
    // auto end = high_res_clock::now();
    // std::chrono::milliseconds diff =
    // std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "radix sort time: " << diff.count() << " ms" << std::endl;
  }
  delete[] control_array;
  delete[] array;
}
void exponential_generator_int64_(int exp_cutoff, double exp_lambda, int line) {
  int n = 1000000000;
  std::vector<std::pair<uintT, uintT>> A(n);
  // pbbs::sequence<uintT> nums(exp_cutoff);
  std::vector<uintT>nums(exp_cutoff);

  std::vector<std::pair<uintT, uintT>> B(n);

  /* 1. making nums[] array */
  parallel_for(int i = 0; i < exp_cutoff; i++)
    nums[i] = (double)n * (exp(-exp_lambda * i) * (1 - exp(-exp_lambda)));

  uintT offset = 0;
  for (uintT i = 0; i < exp_cutoff; i++) {
    offset += nums[i];
  }
  // uintT offset = pbbs::reduce(nums, pbbs::addm<uintT>()); // cout << "offset
  // = " << offset << endl;
  nums[0] += (n - offset);

  uintT offset_ = 0;
  for (uintT i = 0; i < exp_cutoff; i++) {
    offset_ += nums[i];
  }
  // checking if the sum of nums[] equals to n
  if ((uintT)offset_ == (uintT)n) {
    cout << "sum of nums[] == n" << endl;
  }


  /* 2. scan to calculate position */
  uintT *addr = new uintT[exp_cutoff];
  parallel_for (uintT i = 0; i < exp_cutoff; i++)
    addr[i] = nums[i];
  scan_inplace__(addr, exp_cutoff); // store all addresses into addr[]

  /* 3. distribute random numbers into A[i].first */
  parallel_for (uintT i = 0; i < exp_cutoff;i++){
    uintT st = (i == 0) ? 0 : addr[i - 1],
               ed = (i == (uintT)exp_cutoff - 1) ? n : addr[i];
    for (uintT j = st; j < ed; j++) {
          A[j].first = hash64_2(i);
        }
  }
  parallel_for(uintT i = 0; i < n; i++)
    A[i].second = hash64_2(i);

  /* 4. shuffle the keys */
  random_shuffle(A.begin(), A.end());

  delete[] addr;

  uintT length = n;
  uintTPair *control_array = newA(uintTPair, length);
  parallel_for(uintT i = 0; i < length; i++) {
    control_array[i].first = A[i].first;
    control_array[i].second = A[i].second;
  }
  uintTPair *array = newA(uintTPair, length);
  std::cout << "runing test of exponetial input with line " << line << std::endl;
  for (uintT round = 0; round < 5; round++) {
    parallel_for(uintT i = 0; i < length; i++) {
      array[i].first = control_array[i].first;
      array[i].second = control_array[i].second;
    }
    timer t;
    t.start();
    // auto start = high_res_clock::now();
    parallelIntegerSort(array, length, utils::firstF<uintT, uintT>());
    t.stop();
    cout << t.get_total() << endl;
    check_sorted_pairs(array, length, utils::firstF<uintT, uintT>());
    // check_sorted_pairs(array, length, utils::firstF<uintT, uintT>());
    // auto end = high_res_clock::now();
    // std::chrono::milliseconds diff =
    // std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "radix sort time: " << diff.count() << " ms" << std::endl;
  }
  delete[] control_array;
  delete[] array;
}
void test_uniform(uint64_t uniform_max_range, int line) {
  int n = 1000000000;
  std::vector<std::pair<uintT, uintT>> A(n);
  parallel_for(uintT i = 0; i < n; i++) {
    A[i].first = hash64_2(i) % uniform_max_range;
    if (A[i].first > uniform_max_range)
      A[i].first -= uniform_max_range;
    if (A[i].first > uniform_max_range)
      std::cout << "wrong..." << std::endl;
    A[i].first = hash64_2(A[i].first);
    A[i].second = hash64_2(i);
  }
  uintT length = n;
  uintTPair *control_array = newA(uintTPair, length);
  parallel_for(uintT i = 0; i < length; i++) {
    control_array[i].first = A[i].first;
    control_array[i].second = A[i].second;
  }
  uintTPair *array = newA(uintTPair, length);
  std::cout << "runing test of uniform input with line " << line << std::endl;
  for (uintT round = 0; round < 5; round++) {
    parallel_for(uintT i = 0; i < length; i++) {
      array[i].first = control_array[i].first;
      array[i].second = control_array[i].second;
    }
    timer t;
    t.start();
    // auto start = high_res_clock::now();
    parallelIntegerSort(array, length, utils::firstF<uintT, uintT>());
    t.stop();
    cout << t.get_total() << endl;
    // check_sorted_pairs(array, length, utils::firstF<uintT, uintT>());
    // auto end = high_res_clock::now();
    // std::chrono::milliseconds diff =
    // std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "radix sort time: " << diff.count() << " ms" << std::endl;
  }
  delete[] control_array;
  delete[] array;
}

void testgraph(std::string graph_name) {
  std::ifstream fin;
  fin.open(graph_name);
  uintT length;
  fin >> length;
  uintTPair *control_array = newA(uintTPair, length);
  for (uintT i = 0; i < length; i++)
    fin >> control_array[i].first >> control_array[i].second;
  std::cout << "running test with input file " << graph_name << std::endl;
  uintTPair *array = newA(uintTPair, length);
  for (uintT round = 0; round < 5; round++) {
    parallel_for(uintT i = 0; i < length; i++) {
      array[i].first = control_array[i].first;
      array[i].second = control_array[i].second;
    }
    timer t;
    t.start();
    // auto start = high_res_clock::now();
    parallelIntegerSort(array, length, utils::firstF<uintT, uintT>());
    t.stop();
    cout << t.get_total() << endl;
    check_sorted_pairs(array, length, utils::firstF<uintT, uintT>());
    // auto end = high_res_clock::now();
    // std::chrono::milliseconds diff =
    //     std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "radix sort time: " << diff.count() << " ms" << std::endl;
  }
  delete[] control_array;
  delete[] array;
}

int main(int argc, char **argv) {
  std::cout << "the size of uintT is " << sizeof(uintT) << std::endl;
  int id = atoi(argv[1]);
  if (id > 3 && id < 19)
    test_uniform(uniform[id], id);
  else if (id >20 && id < 28) {
    exponential_generator_int64_(1000000, exp_lambda[id-21], id);
  }
  else if (id > 28 && id < 36)
    zipfian_generator_int64_(zipfan[id - 29], id);
  else {
    switch (id) {
    case 37:
      testgraph("/data/zwang358/semisort/graph/com-orkut.txt");
      break;
    case 38:
      testgraph("/data/zwang358/semisort/graph/twitter.txt");
      break;
    case 39:
      testgraph("/data/zwang358/semisort/graph/yahoo_g9.txt");
      break;
    case 40:
      testgraph("/data/zwang358/semisort/graph/sd_arc.txt");
      break;
    case 41:
      testgraph("/data/zwang358/semisort/graph/enwiki.txt");
      break;
    case 42:
      testgraph("/data/zwang358/semisort/graph/webbase2001.txt");
      break;
    case 43:
      testgraph("/data/zwang358/semisort/graph/uk2002.txt");
      break;
    default:
      break;
    }
  }
  return 0;
}

// int main(int argc, char **argv) {
//   commandLine P(argc, argv, "[-r <rounds>] [-c] <inFile>");
//   char *iFile = P.getArgument(0);
//   int rounds = P.getOptionIntValue("-r", 3);
//   int workers = P.getOptionIntValue("-w", 0);
//   if (workers != 0) {
//     setWorkers(workers);
//   }
//   bool check = P.getOption("-c");

//   cout << "workers = " << getWorkers() << endl;
//   seqData D = readSequenceFromFile(iFile);
//   elementType dt = D.dt;
//   long length = D.n;

//   if (dt == intType) {
//     uintT *array = (uintT *)D.A;
//     uintT *control_array = newA(uintT, length);

//     parallel_for(long i = 0; i < length; i++) control_array[i] = array[i];
//     for (long round = 0; round < rounds; round++) {
//       parallel_for(long i = 0; i < length; i++) array[i] = control_array[i];
//       auto start = high_res_clock::now();
//       parallelIntegerSort(array, length, utils::identityF<uintT>());
//       auto end = high_res_clock::now();
//       std::chrono::milliseconds diff =
//           std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//       std::cout << "radix sort time: " << diff.count() << " ms" << std::endl;
//     }

//     if (check) {
//       check_sorted(array, length, utils::identityF<uintT>());
//     }
//     free(array);
//     free(control_array);
//   } else if (dt == intPairT) {

//     uintTPair *array = (uintTPair *)D.A;

//     uintTPair *control_array = newA(uintTPair, length);
//     parallel_for(long i = 0; i < length; i++) control_array[i] = array[i];

//     for (long round = 0; round < rounds; round++) {
//       parallel_for(long i = 0; i < length; i++) array[i] = control_array[i];
//       auto start = high_res_clock::now();
//       parallelIntegerSort(array, length, utils::firstF<uintT, uintT>());
//       auto end = high_res_clock::now();
//       std::chrono::milliseconds diff =
//           std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
//       std::cout << "radix sort time: " << diff.count() << " ms" << std::endl;
//     }

//     if (check) {
//       check_sorted_pairs(array, length, utils::firstF<uintT, uintT>());
//     };
//     free(array);
//     free(control_array);
//   } else if (dt == intStringPairT) {
//     cout << "to implement" << endl;
//   } else {
//     cout << "input file not of right type" << endl;
//   }

//   return 0;
// }
