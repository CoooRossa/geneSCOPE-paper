// 9.ModuleQualityPerm.cpp (2026-01-08)
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <random>
using namespace Rcpp;

struct WelfordState {
    std::vector<double> mean;
    std::vector<double> m2;
    std::vector<int> n;
};

static inline void welford_update(WelfordState &state, int idx, double value) {
    int n1 = state.n[idx] + 1;
    double delta = value - state.mean[idx];
    double mean = state.mean[idx] + delta / n1;
    double m2 = state.m2[idx] + delta * (value - mean);
    state.n[idx] = n1;
    state.mean[idx] = mean;
    state.m2[idx] = m2;
}

static inline void welford_finalize(const WelfordState &state,
                                    NumericVector &mean_out,
                                    NumericVector &sd_out) {
    int n = state.mean.size();
    mean_out = NumericVector(n, NA_REAL);
    sd_out = NumericVector(n, NA_REAL);
    for (int i = 0; i < n; ++i) {
        if (state.n[i] == 0) continue;
        mean_out[i] = state.mean[i];
        if (state.n[i] > 1) {
            sd_out[i] = std::sqrt(state.m2[i] / (state.n[i] - 1));
        }
    }
}

// [[Rcpp::export]]
Rcpp::List module_quality_perm_null(const IntegerVector membership,
                                    const IntegerVector u,
                                    const IntegerVector v,
                                    const NumericVector w,
                                    const int n_perm = 50,
                                    const int seed = 1,
                                    const bool do_member = true,
                                    const bool do_module = true,
                                    const bool do_conductance = true) {
    const int n_nodes = membership.size();
    if (n_nodes == 0 || n_perm <= 0) {
        return Rcpp::List::create(
            _["member_margin_mean"] = NumericVector(),
            _["member_margin_sd"] = NumericVector(),
            _["module_sep_mean"] = NumericVector(),
            _["module_sep_sd"] = NumericVector(),
            _["module_cond_mean"] = NumericVector(),
            _["module_cond_sd"] = NumericVector()
        );
    }

    int max_label = 0;
    for (int i = 0; i < n_nodes; ++i) {
        if (membership[i] > max_label) max_label = membership[i];
    }
    const int n_modules = max_label;

    std::vector<int> labels(n_nodes);
    for (int i = 0; i < n_nodes; ++i) {
        labels[i] = membership[i] - 1;
    }

    WelfordState member_state;
    WelfordState module_sep_state;
    WelfordState module_cond_state;

    if (do_member) {
        member_state.mean.assign(n_nodes, 0.0);
        member_state.m2.assign(n_nodes, 0.0);
        member_state.n.assign(n_nodes, 0);
    }
    if (do_module) {
        module_sep_state.mean.assign(n_modules, 0.0);
        module_sep_state.m2.assign(n_modules, 0.0);
        module_sep_state.n.assign(n_modules, 0);
        if (do_conductance) {
            module_cond_state.mean.assign(n_modules, 0.0);
            module_cond_state.m2.assign(n_modules, 0.0);
            module_cond_state.n.assign(n_modules, 0);
        }
    }

    std::mt19937 rng(seed);
    const int n_edges = u.size();

    for (int b = 0; b < n_perm; ++b) {
        std::vector<int> perm_labels = labels;
        for (int i = n_nodes - 1; i > 0; --i) {
            std::uniform_int_distribution<int> dist(0, i);
            int j = dist(rng);
            std::swap(perm_labels[i], perm_labels[j]);
        }

        std::vector<double> within_sum_node(n_nodes, 0.0);
        std::vector<int> within_n_node(n_nodes, 0);
        std::vector<double> between_sum_node(n_nodes, 0.0);
        std::vector<int> between_n_node(n_nodes, 0);

        std::vector<double> within_sum_mod(n_modules, 0.0);
        std::vector<int> within_n_mod(n_modules, 0);
        std::vector<double> between_sum_mod(n_modules, 0.0);
        std::vector<int> between_n_mod(n_modules, 0);

        for (int e = 0; e < n_edges; ++e) {
            int a = u[e];
            int bnode = v[e];
            if (a < 0 || a >= n_nodes || bnode < 0 || bnode >= n_nodes) continue;
            double wt = w[e];
            if (!R_finite(wt)) continue;
            int ma = perm_labels[a];
            int mb = perm_labels[bnode];
            if (ma == mb) {
                within_sum_node[a] += wt;
                within_sum_node[bnode] += wt;
                within_n_node[a] += 1;
                within_n_node[bnode] += 1;
                within_sum_mod[ma] += wt;
                within_n_mod[ma] += 1;
            } else {
                between_sum_node[a] += wt;
                between_sum_node[bnode] += wt;
                between_n_node[a] += 1;
                between_n_node[bnode] += 1;
                between_sum_mod[ma] += wt;
                between_sum_mod[mb] += wt;
                between_n_mod[ma] += 1;
                between_n_mod[mb] += 1;
            }
        }

        if (do_member) {
            for (int i = 0; i < n_nodes; ++i) {
                if (within_n_node[i] == 0 || between_n_node[i] == 0) continue;
                double within_mean = within_sum_node[i] / within_n_node[i];
                double between_mean = between_sum_node[i] / between_n_node[i];
                double margin = within_mean - between_mean;
                if (R_finite(margin)) {
                    welford_update(member_state, i, margin);
                }
            }
        }

        if (do_module) {
            for (int k = 0; k < n_modules; ++k) {
                if (within_n_mod[k] == 0 || between_n_mod[k] == 0) continue;
                double cohesion = within_sum_mod[k] / within_n_mod[k];
                double between_mean = between_sum_mod[k] / between_n_mod[k];
                double separation = cohesion - between_mean;
                if (R_finite(separation)) {
                    welford_update(module_sep_state, k, separation);
                }
                if (do_conductance) {
                    double denom = within_sum_mod[k] + between_sum_mod[k];
                    if (denom > 0) {
                        double cond = between_sum_mod[k] / denom;
                        if (R_finite(cond)) {
                            welford_update(module_cond_state, k, cond);
                        }
                    }
                }
            }
        }
    }

    NumericVector member_margin_mean;
    NumericVector member_margin_sd;
    NumericVector module_sep_mean;
    NumericVector module_sep_sd;
    NumericVector module_cond_mean;
    NumericVector module_cond_sd;

    if (do_member) {
        welford_finalize(member_state, member_margin_mean, member_margin_sd);
    }
    if (do_module) {
        welford_finalize(module_sep_state, module_sep_mean, module_sep_sd);
        if (do_conductance) {
            welford_finalize(module_cond_state, module_cond_mean, module_cond_sd);
        }
    }

    return Rcpp::List::create(
        _["member_margin_mean"] = member_margin_mean,
        _["member_margin_sd"] = member_margin_sd,
        _["module_sep_mean"] = module_sep_mean,
        _["module_sep_sd"] = module_sep_sd,
        _["module_cond_mean"] = module_cond_mean,
        _["module_cond_sd"] = module_cond_sd
    );
}
