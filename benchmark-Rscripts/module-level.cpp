// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#if defined(__APPLE__) && !defined(RCPP_PARALLEL_USE_TBB)
#define RCPP_PARALLEL_USE_TBB 0
#endif
#include <RcppParallel.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <limits>

using namespace Rcpp;
using namespace RcppParallel;

struct EdgeMap {
  std::unordered_map<uint64_t, double> map;
};

inline uint64_t pack_pair(uint32_t a, uint32_t b) {
  return (static_cast<uint64_t>(a) << 32) | static_cast<uint64_t>(b);
}

inline uint64_t mix_seed(uint64_t seed, uint64_t i) {
  uint64_t x = seed + 0x9E3779B97F4A7C15ULL + i * 0xBF58476D1CE4E5B9ULL;
  x ^= x >> 30;
  x *= 0xBF58476D1CE4E5B9ULL;
  x ^= x >> 27;
  x *= 0x94D049BB133111EBULL;
  x ^= x >> 31;
  return x;
}

class XorShift64 {
 public:
  explicit XorShift64(uint64_t seed) {
    state_ = seed ? seed : 0x106689D45497FDB5ULL;
  }

  inline uint64_t next() {
    uint64_t x = state_;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    state_ = x;
    return x * 2685821657736338717ULL;
  }

  inline uint64_t uniform_u64(uint64_t bound) {
    if (bound == 0) return 0;
    uint64_t x = next();
    uint64_t lim = std::numeric_limits<uint64_t>::max();
    uint64_t m = lim - (lim % bound);
    while (x >= m) x = next();
    return x % bound;
  }

  inline int uniform_int(int low, int high) {
    uint64_t bound = static_cast<uint64_t>(high - low + 1);
    return low + static_cast<int>(uniform_u64(bound));
  }

 private:
  uint64_t state_;
};

struct Sampler {
  std::vector<int> marks;
  int mark_id;
  std::vector<int> selected_positions;

  explicit Sampler(std::size_t n) : marks(n, 0), mark_id(1) {}

  void sample(int k, XorShift64 &rng, const std::vector<int> &bg, std::vector<int> &out) {
    out.clear();
    if (k <= 0) return;
    int n = static_cast<int>(bg.size());
    if (k > n) return;

    selected_positions.clear();
    selected_positions.reserve(k);

    for (int j = n - k; j < n; ++j) {
      int t = rng.uniform_int(0, j);
      int pick = (marks[t] == mark_id) ? j : t;
      marks[pick] = mark_id;
      selected_positions.push_back(pick);
    }

    out.reserve(k);
    for (int pos : selected_positions) {
      out.push_back(bg[pos]);
    }
    std::sort(out.begin(), out.end());

    mark_id++;
    if (mark_id == std::numeric_limits<int>::max()) {
      std::fill(marks.begin(), marks.end(), 0);
      mark_id = 1;
    }
  }
};

inline void sum_pairs_sorted(const std::vector<int> &genes, EdgeMap *edge_map,
                             int &n_pairs_mapped, double &sum_mapped) {
  n_pairs_mapped = 0;
  sum_mapped = 0.0;
  int k = static_cast<int>(genes.size());
  if (k < 2) return;

  for (int a = 0; a < k; ++a) {
    uint32_t ga = static_cast<uint32_t>(genes[a]);
    for (int b = a + 1; b < k; ++b) {
      uint32_t gb = static_cast<uint32_t>(genes[b]);
      uint64_t key = pack_pair(ga, gb);
      auto it = edge_map->map.find(key);
      if (it != edge_map->map.end()) {
        n_pairs_mapped++;
        sum_mapped += it->second;
      }
    }
  }
}

inline void pair_stats_sorted_threshold(const std::vector<int> &genes, EdgeMap *edge_map,
                                        double threshold, int &n_pairs_mapped, double &sum_mapped,
                                        int &n_pairs_ge_threshold) {
  n_pairs_mapped = 0;
  sum_mapped = 0.0;
  n_pairs_ge_threshold = 0;
  int k = static_cast<int>(genes.size());
  if (k < 2) return;

  for (int a = 0; a < k; ++a) {
    uint32_t ga = static_cast<uint32_t>(genes[a]);
    for (int b = a + 1; b < k; ++b) {
      uint32_t gb = static_cast<uint32_t>(genes[b]);
      uint64_t key = pack_pair(ga, gb);
      auto it = edge_map->map.find(key);
      if (it == edge_map->map.end()) continue;
      n_pairs_mapped++;
      double score = it->second;
      sum_mapped += score;
      if (score >= threshold) n_pairs_ge_threshold++;
    }
  }
}

// [[Rcpp::export]]
SEXP build_edge_map_cpp(IntegerVector edge_i, IntegerVector edge_j, NumericVector score) {
  XPtr<EdgeMap> ptr(new EdgeMap(), true);
  std::size_t n = static_cast<std::size_t>(edge_i.size());
  ptr->map.reserve(static_cast<std::size_t>(n * 1.3) + 1);

  for (int idx = 0; idx < edge_i.size(); ++idx) {
    int i = edge_i[idx];
    int j = edge_j[idx];
    if (i < 0 || j < 0) continue;
    if (i == j) continue;
    uint32_t a = static_cast<uint32_t>(i);
    uint32_t b = static_cast<uint32_t>(j);
    if (a > b) std::swap(a, b);
    uint64_t key = pack_pair(a, b);
    double s = score[idx];
    auto it = ptr->map.find(key);
    if (it == ptr->map.end() || s > it->second) {
      ptr->map[key] = s;
    }
  }

  return ptr;
}

// [[Rcpp::export]]
List module_pair_stats_cpp(IntegerVector genes, SEXP edge_map_ptr) {
  XPtr<EdgeMap> edge_map(edge_map_ptr);
  std::vector<int> gene_vec;
  gene_vec.reserve(genes.size());
  for (int g : genes) {
    if (g >= 0) gene_vec.push_back(g);
  }
  if (gene_vec.size() < 2) {
    return List::create(
      _["n_pairs_mapped"] = 0,
      _["sum_mapped"] = 0.0
    );
  }
  std::sort(gene_vec.begin(), gene_vec.end());
  gene_vec.erase(std::unique(gene_vec.begin(), gene_vec.end()), gene_vec.end());

  int n_pairs_mapped = 0;
  double sum_mapped = 0.0;
  sum_pairs_sorted(gene_vec, edge_map.get(), n_pairs_mapped, sum_mapped);

  return List::create(
    _["n_pairs_mapped"] = n_pairs_mapped,
    _["sum_mapped"] = sum_mapped
  );
}

// [[Rcpp::export]]
List module_pair_stats_threshold_cpp(IntegerVector genes, SEXP edge_map_ptr, double threshold) {
  XPtr<EdgeMap> edge_map(edge_map_ptr);
  std::vector<int> gene_vec;
  gene_vec.reserve(genes.size());
  for (int g : genes) {
    if (g >= 0) gene_vec.push_back(g);
  }
  if (gene_vec.size() < 2) {
    return List::create(
      _["n_pairs_mapped"] = 0,
      _["sum_mapped"] = 0.0,
      _["n_pairs_ge_threshold"] = 0
    );
  }
  std::sort(gene_vec.begin(), gene_vec.end());
  gene_vec.erase(std::unique(gene_vec.begin(), gene_vec.end()), gene_vec.end());

  int n_pairs_mapped = 0;
  double sum_mapped = 0.0;
  int n_pairs_ge_threshold = 0;
  pair_stats_sorted_threshold(gene_vec, edge_map.get(), threshold,
                              n_pairs_mapped, sum_mapped, n_pairs_ge_threshold);

  return List::create(
    _["n_pairs_mapped"] = n_pairs_mapped,
    _["sum_mapped"] = sum_mapped,
    _["n_pairs_ge_threshold"] = n_pairs_ge_threshold
  );
}

struct NullEavgWorker : public Worker {
  const std::vector<int> &bg;
  const std::vector<int> &module;
  int module_size;
  double denom;
  EdgeMap *edge_map;
  uint64_t seed;
  RVector<double> out;

  NullEavgWorker(const std::vector<int> &bg, const std::vector<int> &module, int module_size,
                 double denom, EdgeMap *edge_map, uint64_t seed, NumericVector out)
    : bg(bg),
      module(module),
      module_size(module_size),
      denom(denom),
      edge_map(edge_map),
      seed(seed),
      out(out) {}

  void operator()(std::size_t begin, std::size_t end) {
    Sampler sampler(bg.size());
    std::vector<int> sampled;
    std::vector<int> genes;
    for (std::size_t i = begin; i < end; ++i) {
      double eavg = NA_REAL;
      if (module_size < 2) {
        out[i] = eavg;
        continue;
      }
      uint64_t local_seed = mix_seed(seed, static_cast<uint64_t>(i + 1));
      XorShift64 rng(local_seed);
      int anchor_idx = rng.uniform_int(0, module_size - 1);
      int anchor = module[anchor_idx];
      sampler.sample(module_size - 1, rng, bg, sampled);
      if (sampled.size() != static_cast<std::size_t>(module_size - 1)) {
        out[i] = eavg;
        continue;
      }
      genes = sampled;
      genes.push_back(anchor);
      std::sort(genes.begin(), genes.end());
      int n_pairs_mapped = 0;
      double sum_mapped = 0.0;
      sum_pairs_sorted(genes, edge_map, n_pairs_mapped, sum_mapped);
      eavg = (denom > 0.0) ? (sum_mapped / denom) : NA_REAL;
      out[i] = eavg;
    }
  }
};

// [[Rcpp::export]]
NumericVector draw_null_eavg_cpp(int n_random, IntegerVector module_indices,
                                 IntegerVector bg_indices, SEXP edge_map_ptr, double seed) {
  if (n_random <= 0) return NumericVector(0);
  NumericVector out(n_random, NA_REAL);
  std::vector<int> module(module_indices.begin(), module_indices.end());
  module.erase(std::remove_if(module.begin(), module.end(), [](int v) { return v < 0; }), module.end());
  std::sort(module.begin(), module.end());
  module.erase(std::unique(module.begin(), module.end()), module.end());
  int module_size = static_cast<int>(module.size());
  if (module_size <= 0) return out;

  std::vector<int> bg(bg_indices.begin(), bg_indices.end());
  bg.erase(std::remove_if(bg.begin(), bg.end(), [](int v) { return v < 0; }), bg.end());
  if (bg.size() < static_cast<std::size_t>(module_size - 1)) return out;

  XPtr<EdgeMap> edge_map(edge_map_ptr);
  double denom = static_cast<double>(module_size) * (module_size - 1) / 2.0;
  NullEavgWorker worker(bg, module, module_size, denom, edge_map.get(),
                        static_cast<uint64_t>(seed), out);
  parallelFor(0, static_cast<std::size_t>(n_random), worker);
  return out;
}

struct NullEavgPosFracWorker : public Worker {
  const std::vector<int> &bg;
  const std::vector<int> &module;
  int module_size;
  double denom;
  double threshold;
  EdgeMap *edge_map;
  uint64_t seed;
  RVector<double> out_eavg;
  RVector<double> out_posfrac;

  NullEavgPosFracWorker(const std::vector<int> &bg, const std::vector<int> &module, int module_size,
                        double denom, double threshold, EdgeMap *edge_map, uint64_t seed,
                        NumericVector out_eavg, NumericVector out_posfrac)
    : bg(bg),
      module(module),
      module_size(module_size),
      denom(denom),
      threshold(threshold),
      edge_map(edge_map),
      seed(seed),
      out_eavg(out_eavg),
      out_posfrac(out_posfrac) {}

  void operator()(std::size_t begin, std::size_t end) {
    Sampler sampler(bg.size());
    std::vector<int> sampled;
    std::vector<int> genes;
    for (std::size_t i = begin; i < end; ++i) {
      double eavg = NA_REAL;
      double posfrac = NA_REAL;
      if (module_size < 2) {
        out_eavg[i] = eavg;
        out_posfrac[i] = posfrac;
        continue;
      }
      uint64_t local_seed = mix_seed(seed, static_cast<uint64_t>(i + 1));
      XorShift64 rng(local_seed);
      int anchor_idx = rng.uniform_int(0, module_size - 1);
      int anchor = module[anchor_idx];
      sampler.sample(module_size - 1, rng, bg, sampled);
      if (sampled.size() != static_cast<std::size_t>(module_size - 1)) {
        out_eavg[i] = eavg;
        out_posfrac[i] = posfrac;
        continue;
      }
      genes = sampled;
      genes.push_back(anchor);
      std::sort(genes.begin(), genes.end());
      int n_pairs_mapped = 0;
      double sum_mapped = 0.0;
      int n_pairs_ge_threshold = 0;
      pair_stats_sorted_threshold(genes, edge_map, threshold,
                                  n_pairs_mapped, sum_mapped, n_pairs_ge_threshold);
      if (denom > 0.0) {
        eavg = sum_mapped / denom;
        posfrac = static_cast<double>(n_pairs_ge_threshold) / denom;
      }
      out_eavg[i] = eavg;
      out_posfrac[i] = posfrac;
    }
  }
};

// [[Rcpp::export]]
List draw_null_eavg_posfrac_cpp(int n_random, IntegerVector module_indices,
                                IntegerVector bg_indices, SEXP edge_map_ptr,
                                double threshold, double seed) {
  if (n_random <= 0) {
    return List::create(
      _["null_eavg"] = NumericVector(0),
      _["null_pos_pair_frac"] = NumericVector(0)
    );
  }
  NumericVector out_eavg(n_random, NA_REAL);
  NumericVector out_posfrac(n_random, NA_REAL);
  std::vector<int> module(module_indices.begin(), module_indices.end());
  module.erase(std::remove_if(module.begin(), module.end(), [](int v) { return v < 0; }), module.end());
  std::sort(module.begin(), module.end());
  module.erase(std::unique(module.begin(), module.end()), module.end());
  int module_size = static_cast<int>(module.size());
  if (module_size <= 0) {
    return List::create(
      _["null_eavg"] = out_eavg,
      _["null_pos_pair_frac"] = out_posfrac
    );
  }

  std::vector<int> bg(bg_indices.begin(), bg_indices.end());
  bg.erase(std::remove_if(bg.begin(), bg.end(), [](int v) { return v < 0; }), bg.end());
  if (bg.size() < static_cast<std::size_t>(module_size - 1)) {
    return List::create(
      _["null_eavg"] = out_eavg,
      _["null_pos_pair_frac"] = out_posfrac
    );
  }

  XPtr<EdgeMap> edge_map(edge_map_ptr);
  double denom = static_cast<double>(module_size) * (module_size - 1) / 2.0;
  NullEavgPosFracWorker worker(bg, module, module_size, denom, threshold, edge_map.get(),
                               static_cast<uint64_t>(seed), out_eavg, out_posfrac);
  parallelFor(0, static_cast<std::size_t>(n_random), worker);
  return List::create(
    _["null_eavg"] = out_eavg,
    _["null_pos_pair_frac"] = out_posfrac
  );
}

struct NullAnchorCrossEavgPosFracWorker : public Worker {
  const std::vector<int> &module;
  const std::vector<int> &bg;
  int module_size;
  int bg_size;
  int n_edges_target;
  double denom;
  double threshold;
  EdgeMap *edge_map;
  uint64_t seed;
  RVector<double> out_eavg;
  RVector<double> out_posfrac;

  NullAnchorCrossEavgPosFracWorker(const std::vector<int> &module, const std::vector<int> &bg,
                                  int module_size, int bg_size, int n_edges_target,
                                  double denom, double threshold, EdgeMap *edge_map,
                                  uint64_t seed, NumericVector out_eavg, NumericVector out_posfrac)
    : module(module),
      bg(bg),
      module_size(module_size),
      bg_size(bg_size),
      n_edges_target(n_edges_target),
      denom(denom),
      threshold(threshold),
      edge_map(edge_map),
      seed(seed),
      out_eavg(out_eavg),
      out_posfrac(out_posfrac) {}

  void operator()(std::size_t begin, std::size_t end) {
    const uint64_t N = static_cast<uint64_t>(module_size) * static_cast<uint64_t>(bg_size);
    std::unordered_map<uint64_t, uint64_t> swaps;
    swaps.reserve(static_cast<std::size_t>(n_edges_target) * 2 + 16);

    for (std::size_t i = begin; i < end; ++i) {
      if (module_size <= 0 || bg_size <= 0 || n_edges_target <= 0 || denom <= 0.0 || N == 0) {
        out_eavg[i] = NA_REAL;
        out_posfrac[i] = NA_REAL;
        continue;
      }
      uint64_t local_seed = mix_seed(seed, static_cast<uint64_t>(i + 1));
      XorShift64 rng(local_seed);
      swaps.clear();

      double sum = 0.0;
      int pos = 0;

      for (uint64_t draw = 0; draw < static_cast<uint64_t>(n_edges_target); ++draw) {
        uint64_t remaining = N - draw;
        uint64_t j = rng.uniform_u64(remaining);

        uint64_t idx = j;
        auto it = swaps.find(j);
        if (it != swaps.end()) idx = it->second;

        uint64_t last = remaining - 1;
        uint64_t last_val = last;
        auto it_last = swaps.find(last);
        if (it_last != swaps.end()) last_val = it_last->second;
        swaps[j] = last_val;

        uint64_t mod_pos = idx % static_cast<uint64_t>(module_size);
        uint64_t bg_pos = idx / static_cast<uint64_t>(module_size);
        uint32_t u = static_cast<uint32_t>(module[static_cast<std::size_t>(mod_pos)]);
        uint32_t v = static_cast<uint32_t>(bg[static_cast<std::size_t>(bg_pos)]);
        if (u > v) std::swap(u, v);
        uint64_t key = pack_pair(u, v);
        auto it_edge = edge_map->map.find(key);
        double score = (it_edge == edge_map->map.end()) ? 0.0 : it_edge->second;
        sum += score;
        if (score >= threshold) pos++;
      }

      out_eavg[i] = sum / denom;
      out_posfrac[i] = static_cast<double>(pos) / denom;
    }
  }
};

// [[Rcpp::export]]
List draw_null_anchor_cross_eavg_posfrac_cpp(int n_random, IntegerVector module_indices,
                                             IntegerVector bg_indices, int n_edges_target,
                                             SEXP edge_map_ptr, double threshold, int seed) {
  if (n_random <= 0 || n_edges_target <= 0) {
    return List::create(
      _["null_eavg"] = NumericVector(0),
      _["null_pos_pair_frac"] = NumericVector(0)
    );
  }

  std::vector<int> module(module_indices.begin(), module_indices.end());
  module.erase(std::remove_if(module.begin(), module.end(), [](int v) { return v < 0; }), module.end());
  std::sort(module.begin(), module.end());
  module.erase(std::unique(module.begin(), module.end()), module.end());
  int module_size = static_cast<int>(module.size());

  std::vector<int> bg(bg_indices.begin(), bg_indices.end());
  bg.erase(std::remove_if(bg.begin(), bg.end(), [](int v) { return v < 0; }), bg.end());
  std::sort(bg.begin(), bg.end());
  bg.erase(std::unique(bg.begin(), bg.end()), bg.end());
  int bg_size = static_cast<int>(bg.size());

  uint64_t N = static_cast<uint64_t>(module_size) * static_cast<uint64_t>(bg_size);
  if (module_size <= 0 || bg_size <= 0 || N < static_cast<uint64_t>(n_edges_target)) {
    return List::create(
      _["null_eavg"] = NumericVector(0),
      _["null_pos_pair_frac"] = NumericVector(0)
    );
  }

  NumericVector out_eavg(n_random, NA_REAL);
  NumericVector out_posfrac(n_random, NA_REAL);

  XPtr<EdgeMap> edge_map(edge_map_ptr);
  double denom = static_cast<double>(n_edges_target);
  NullAnchorCrossEavgPosFracWorker worker(
    module, bg, module_size, bg_size, n_edges_target, denom, threshold,
    edge_map.get(), static_cast<uint64_t>(seed), out_eavg, out_posfrac
  );
  parallelFor(0, static_cast<std::size_t>(n_random), worker);

  return List::create(
    _["null_eavg"] = out_eavg,
    _["null_pos_pair_frac"] = out_posfrac
  );
}
