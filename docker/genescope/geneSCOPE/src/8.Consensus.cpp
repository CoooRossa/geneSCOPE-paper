// Consensus accelerators: build sparse co-occurrence (COO) from membership runs
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <cstdint>

using namespace Rcpp;

#ifdef _OPENMP
  #include <omp.h>
#endif

// 64-bit key for an unordered pair (i<j) with 0-based indices
static inline uint64_t pair_key32(uint32_t a, uint32_t b) {
  if (a > b) std::swap(a, b);
  return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

// Exported name intentionally suffixed with _cpp so that Rcpp::sourceCpp()
// does not clash with an R wrapper named consensus_coo().
// [[Rcpp::export(name="consensus_coo_cpp")]]
Rcpp::List consensus_coo_cpp(Rcpp::IntegerMatrix memb,
                         double thr = 0.0,
                         int n_threads = 1) {
  const int N = memb.nrow();
  const int R = memb.ncol();
  if (N <= 1 || R <= 0) {
    return Rcpp::List::create(
      Rcpp::Named("i") = Rcpp::IntegerVector(0),
      Rcpp::Named("j") = Rcpp::IntegerVector(0),
      Rcpp::Named("x") = Rcpp::NumericVector(0),
      Rcpp::Named("n") = N,
      Rcpp::Named("runs") = R
    );
  }

  // Minimum co-occurrence count required (strictly positive when thr==0)
  int kmin;
  if (thr <= 0.0) kmin = 1; else kmin = (int)std::ceil(thr * (double)R - 1e-12);
  if (kmin < 1) kmin = 1;

  // Global accumulator of counts across runs
  std::unordered_map<uint64_t, uint16_t> counts;
  counts.reserve(std::min<size_t>((size_t)N * 8u, 1000000u));

  // OpenMP thread control
  #ifdef _OPENMP
    if (n_threads < 1) n_threads = 1;
  #else
    n_threads = 1;
  #endif

  // Process runs one by one; within each run parallelise over clusters
  for (int r = 0; r < R; ++r) {
    // Build cluster -> member indices map for this run
    // Collect unique labels first to avoid many small map operations
    std::unordered_map<int, int> lab2idx;
    lab2idx.reserve((size_t)N / 8 + 8);
    std::vector< std::vector<int> > groups;
    groups.reserve((size_t)N / 8 + 8);

    for (int i = 0; i < N; ++i) {
      int lab = memb(i, r);
      auto it = lab2idx.find(lab);
      if (it == lab2idx.end()) {
        int nid = (int)groups.size();
        lab2idx.emplace(lab, nid);
        groups.emplace_back();
        groups.back().push_back(i);
      } else {
        groups[it->second].push_back(i);
      }
    }

    // Thread-local maps to avoid contention while iterating pairs
    #ifdef _OPENMP
      std::vector< std::unordered_map<uint64_t, uint16_t> > local;
      local.resize((size_t)n_threads);
      for (int t = 0; t < n_threads; ++t) local[t].reserve(1024);

      #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
      for (int g = 0; g < (int)groups.size(); ++g) {
        int tid = 0;
        #ifdef _OPENMP
          tid = omp_get_thread_num();
        #endif
        auto &acc = local[tid];
        const std::vector<int> &mem = groups[g];
        const int s = (int)mem.size();
        if (s <= 1) continue;
        for (int a = 0; a < s - 1; ++a) {
          const uint32_t ia = (uint32_t)mem[a];
          for (int b = a + 1; b < s; ++b) {
            const uint32_t ib = (uint32_t)mem[b];
            const uint64_t key = pair_key32(ia, ib);
            auto it2 = acc.find(key);
            if (it2 == acc.end()) acc.emplace(key, (uint16_t)1);
            else {
              uint16_t v = it2->second;
              if (v < (uint16_t)65535) it2->second = (uint16_t)(v + 1);
            }
          }
        }
      }

      // Merge thread-local maps into global counts
      for (int t = 0; t < n_threads; ++t) {
        for (auto &kv : local[t]) {
          auto it = counts.find(kv.first);
          if (it == counts.end()) counts.emplace(kv.first, kv.second);
          else {
            uint32_t sum = (uint32_t)it->second + (uint32_t)kv.second;
            it->second = (uint16_t)std::min<uint32_t>(sum, 65535u);
          }
        }
      }
    #else
      // No OpenMP: single-thread accumulate directly into global map
      for (size_t g = 0; g < groups.size(); ++g) {
        const std::vector<int> &mem = groups[g];
        const int s = (int)mem.size();
        if (s <= 1) continue;
        for (int a = 0; a < s - 1; ++a) {
          const uint32_t ia = (uint32_t)mem[a];
          for (int b = a + 1; b < s; ++b) {
            const uint32_t ib = (uint32_t)mem[b];
            const uint64_t key = pair_key32(ia, ib);
            auto it2 = counts.find(key);
            if (it2 == counts.end()) counts.emplace(key, (uint16_t)1);
            else {
              uint16_t v = it2->second;
              if (v < (uint16_t)65535) it2->second = (uint16_t)(v + 1);
            }
          }
        }
      }
    #endif
  }

  // Emit COO arrays for pairs meeting threshold
  std::vector<int> I; I.reserve(counts.size());
  std::vector<int> J; J.reserve(counts.size());
  std::vector<double> X; X.reserve(counts.size());
  for (const auto &kv : counts) {
    const uint16_t c = kv.second;
    if ((int)c >= kmin) {
      const uint32_t a = (uint32_t)(kv.first >> 32);
      const uint32_t b = (uint32_t)(kv.first & 0xffffffffu);
      I.push_back((int)a + 1);
      J.push_back((int)b + 1);
      X.push_back((double)c / (double)R);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("i") = Rcpp::wrap(I),
    Rcpp::Named("j") = Rcpp::wrap(J),
    Rcpp::Named("x") = Rcpp::wrap(X),
    Rcpp::Named("n") = N,
    Rcpp::Named("runs") = R
  );
}
