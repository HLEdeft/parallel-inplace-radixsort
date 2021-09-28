#include "get_time.h"
#include "parse_command_line.h"
#include "utilities.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <climits>
#include "sequence.h"
#include "random_shuffle.h"
using namespace std;
using namespace pbbs;
void uniform_generator_int64_(uint64_t uniform_max_range, int line) {
  int n = 1000000000;
  vector<std::pair<uint64_t, uint64_t>> A(n);
  parallel_for(
      0, n,
      [&](uint32_t i) {
        // in order to put all keys in range [0, uniform_max_range]
        A[i].first = pbbs::hash64_2(i) % uniform_max_range;
        if (A[i].first > uniform_max_range)
          A[i].first -= uniform_max_range;
        if (A[i].first > uniform_max_range)
          cout << "wrong..." << endl;
        A[i].first = pbbs::hash64_2(A[i].first);
        A[i].second = pbbs::hash64_2(i);
      },
      1);
  ofstream fout;
  string s = "/data/zwang358/semisort/uniform/" + to_string(line) + ".txt";
  fout.open(s);
  fout << n << endl;
  for (int i = 0; i < n; i++)
    fout << A[i].first << " " << A[i].second << endl;
  fout.close();
}
void scan_inplace__ (uint32_t* in, uint32_t n) {
    if (n <= 1024) {
        // cout << "THRESHOLDS: " << THRESHOLDS << endl;
        for (size_t i = 1;i < n;i++) {
            in[i] += in[i-1];
        }
        return;
    }
    uint32_t root_n = (uint32_t)sqrt(n);// split the array into root n blocks
    uint32_t* offset = new uint32_t[root_n-1];

    parallel_for (0, root_n-1, [&] (size_t i) {
        offset[i] = 0;
        for (size_t j = i*root_n;j < (i+1)*root_n;j++) {
            offset[i] += in[j];
        }
    });

    for (size_t i = 1;i < root_n-1;i++) offset[i] += offset[i-1];

    // prefix sum for each subarray
    parallel_for (0, root_n, [&] (size_t i) {
        if (i == root_n-1) {// the last one
            for (size_t j = i*root_n+1;j < n;j++) {
                in[j] += in[j-1];
            }
        } else {
            for (size_t j = i*root_n+1;j < (i+1)*root_n;j++) {
                in[j] += in[j-1];
            }
        }
    });

    parallel_for (1, root_n, [&] (size_t i) {
        if (i == root_n-1)  {
            for (size_t j = i * root_n;j < n;j++){
                in[j] += offset[i-1];
            }
        } else {
            for (size_t j = i * root_n;j < (i+1)*root_n;j++) {
                in[j] +=  offset[i-1];
            }
        }
    });
    delete[] offset;
}

void zipfian_generator_int64_ (uint32_t zipf_s, int line) {
    int n = 1000000000;
    pbbs::sequence<uint32_t> nums(zipf_s); // in total zipf_s kinds of keys
    pbbs::sequence<std::pair<uint64_t, uint64_t>> B(n);

    /* 1. making nums[] array */
    uint32_t number = (uint32_t) (n / log(n)); // number= n/ln(n)
    parallel_for (0, zipf_s, [&] (uint32_t i) {
        nums[i] = (uint32_t) (number / (i+1));
    }, 1);

    uint32_t offset = pbbs::reduce(nums, pbbs::addm<uint32_t>()); // cout << "offset = " << offset << endl;
    // nums[zipf_s-1] += (n - offset);
    nums[0] += (n - offset);

    // checking if the sum of nums[] equals to n
    if (pbbs::reduce(nums, pbbs::addm<uint32_t>()) == (uint32_t)n) {
        cout << "sum of nums[] == n" << endl;
    }

    /* 2. scan to calculate position */
    // pbbs::sequence<uint32_t> addr(zipf_s);
    uint32_t* addr = new uint32_t[zipf_s];
    parallel_for (0, zipf_s, [&] (uint32_t i) {
        addr[i] = nums[i];
    }, 1);
    scan_inplace__(addr, zipf_s); // store all addresses into addr[]

    /* 3. distribute random numbers into A[i].first */
    parallel_for (0, zipf_s, [&] (size_t i) {
        uint32_t st = (i == 0) ? 0 : addr[i-1],
                 ed = (i == zipf_s-1) ? n : addr[i];
        for (uint32_t j = st; j < ed; j++) {
            B[j].first = pbbs::hash64_2(i);
        }
    }, 1);
    parallel_for (0, n, [&] (size_t i){
        B[i].second = pbbs::hash64_2(i);
    }, 1);

    /* 4. shuffle the keys */
    // random_shuffle(A, A + n);
    pbbs::sequence<std::pair<uint64_t, uint64_t>> C = pbbs::random_shuffle(B, n);

    ofstream fout;
    string s = "/data/zwang358/semisort/zip/" + to_string(line) + ".txt";
    fout.open(s);
    fout << n << endl;
    for (int i = 0; i < n; i++)
      fout << C[i].first << " " << C[i].second << endl;
    fout.close();


    delete[] addr;
}

int main(int argc, char *argv[]) {
  int line = atoi(argv[1]);
  int n = atoi(argv[2]);
  uniform_generator_int64_(n, line);
  zipfian_generator_int64_(n,line);
  return 0;
}



