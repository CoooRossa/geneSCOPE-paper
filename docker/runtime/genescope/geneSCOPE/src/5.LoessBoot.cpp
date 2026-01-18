// LOESS + residual bootstrap with unified extended interface
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD 1
#define ARMA_NO_DEBUG 1
#define ARMA_USE_OPENMP 1
#include <RcppArmadillo.h>
#include <random>
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// --------- Local weighted regression (tricube) ----------
static vec local_loess(const vec &x, const vec &y, const vec &xout,
                       double span, int deg, int k_max)
{
    uword n = x.n_elem, m = xout.n_elem;
    vec fit(m, fill::value(datum::nan));
    if (n == 0)
        return fit;
    if (n < 3)
    {
        fit.fill(mean(y));
        return fit;
    }
    if (span <= 0)
        span = 0.3;
    uword k_raw = std::max<uword>(2, std::ceil(span * n));
    uword k_use = (k_max > 0) ? std::min<uword>(k_raw, (uword)k_max) : k_raw;
    deg = (deg == 0) ? 0 : 1;
    std::vector<uword> order_idx(n);
    std::iota(order_idx.begin(), order_idx.end(), 0);
    for (uword j = 0; j < m; ++j)
    {
        double x0 = xout[j];
        std::nth_element(order_idx.begin(), order_idx.begin() + k_use, order_idx.end(),
                         [&](uword a, uword b)
                         {
                             return std::abs(x[a] - x0) < std::abs(x[b] - x0);
                         });
        vec xs(k_use), ys(k_use);
        for (uword t = 0; t < k_use; ++t)
        {
            xs[t] = x[order_idx[t]];
            ys[t] = y[order_idx[t]];
        }
        double dmax = max(abs(xs - x0));
        if (dmax <= 0)
        {
            fit[j] = mean(ys);
            continue;
        }
        vec w(k_use);
        for (uword t = 0; t < k_use; ++t)
        {
            double u = std::abs((xs[t] - x0) / dmax);
            if (u >= 1)
                w[t] = 0;
            else
            {
                double a = 1 - u * u * u;
                w[t] = a * a * a;
            }
        }
        double W = accu(w);
        if (W <= 0)
        {
            fit[j] = mean(ys);
            continue;
        }
        if (deg == 0)
        {
            fit[j] = dot(w, ys) / W;
        }
        else
        {
            vec Xc = xs - x0;
            double S0 = W;
            double S1 = dot(w, Xc);
            double S2 = dot(w % Xc, Xc);
            double T0 = dot(w, ys);
            double T1 = dot(w % Xc, ys);
            double den = S0 * S2 - S1 * S1;
            if (std::abs(den) < 1e-14)
            {
                fit[j] = T0 / S0;
            }
            else
            {
                double beta1 = (S0 * T1 - S1 * T0) / den;
                double beta0 = (T0 - beta1 * S1) / S0;
                fit[j] = beta0;
            }
        }
    }
    return fit;
}

// Approximate hat diagonal (used only for edf and variance estimates)
static vec hat_diag(const vec &x, double span, int k_max)
{
    uword n = x.n_elem;
    vec h(n, fill::zeros);
    if (n < 3)
        return h;
    if (span <= 0)
        span = 0.3;
    uword k_raw = std::max<uword>(2, std::ceil(span * n));
    uword k_use = (k_max > 0) ? std::min<uword>(k_raw, (uword)k_max) : k_raw;
    std::vector<uword> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    for (uword i = 0; i < n; ++i)
    {
        double x0 = x[i];
        std::nth_element(ord.begin(), ord.begin() + k_use, ord.end(),
                         [&](uword a, uword b)
                         {
                             return std::abs(x[a] - x0) < std::abs(x[b] - x0);
                         });
        double dmax = 0.0;
        for (uword t = 0; t < k_use; ++t)
        {
            double d = std::abs(x[ord[t]] - x0);
            if (d > dmax)
                dmax = d;
        }
        if (dmax <= 0)
        {
            h[i] = 1;
            continue;
        }
        double wsum = 0, selfw = 0;
        for (uword t = 0; t < k_use; ++t)
        {
            uword id = ord[t];
            double u = std::abs((x[id] - x0) / dmax);
            double w = (u >= 1) ? 0 : std::pow(1 - u * u * u, 3);
            if (w > 0)
            {
                wsum += w;
                if (id == i)
                    selfw = w;
            }
        }
        if (wsum > 0)
            h[i] = selfw / wsum;
    }
    return h;
}

// Weight^2 (squared L2 norm of normalized kernel) at grid points
static vec weight_sq(const vec &x, const vec &grid, double span, int k_max)
{
    uword n = x.n_elem, m = grid.n_elem;
    vec out(m, fill::zeros);
    if (n < 3)
        return out;
    if (span <= 0)
        span = 0.3;
    uword k_raw = std::max<uword>(2, std::ceil(span * n));
    uword k_use = (k_max > 0) ? std::min<uword>(k_raw, (uword)k_max) : k_raw;
    std::vector<uword> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    for (uword j = 0; j < m; ++j)
    {
        double x0 = grid[j];
        std::nth_element(ord.begin(), ord.begin() + k_use, ord.end(),
                         [&](uword a, uword b)
                         {
                             return std::abs(x[a] - x0) < std::abs(x[b] - x0);
                         });
        double dmax = 0;
        for (uword t = 0; t < k_use; ++t)
        {
            double d = std::abs(x[ord[t]] - x0);
            if (d > dmax)
                dmax = d;
        }
        if (dmax <= 0)
            continue;
        vec w(k_use);
        for (uword t = 0; t < k_use; ++t)
        {
            double u = std::abs((x[ord[t]] - x0) / dmax);
            if (u >= 1)
                w[t] = 0;
            else
            {
                double a = 1 - u * u * u;
                w[t] = a * a * a;
            }
        }
        double W = accu(w);
        if (W <= 0)
            continue;
        w /= W;
        out[j] = accu(w % w);
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::List loess_residual_bootstrap(const arma::vec &x,
                                    const arma::vec &y,
                                    const arma::uvec &strat,
                                    const arma::vec &grid,
                                    int B = 1000,
                                    double span = 0.45,
                                    int deg = 1,
                                    int n_threads = 1,
                                    int k_max = -1,
                                    bool keep_boot = true,
                                    int adjust_mode = 0,
                                    int ci_type = 0,
                                    double level = 0.95)
{

    if (x.n_elem != y.n_elem || x.n_elem != strat.n_elem)
        stop("Input length mismatch");
    if (x.n_elem < 5)
        stop("Too few points");
    if (level <= 0 || level >= 1)
        stop("level must be in (0,1)");

#ifdef _OPENMP
    if (n_threads < 1)
        n_threads = 1;
    omp_set_num_threads(n_threads);
#endif

    vec fit_grid = local_loess(x, y, grid, span, deg, k_max);
    vec fit_obs = local_loess(x, y, x, span, deg, k_max);
    vec resid = y - fit_obs;

    // strata indexing (1-based expected)
    uword K = strat.max();
    if (K == 0)
        stop("Strata must be 1-based");
    std::vector<uvec> strata(K);
    for (uword k = 1; k <= K; ++k)
    {
        uvec id = find(strat == k);
        strata[k - 1] = id;
    }

    vec h = hat_diag(x, span, k_max);
    double edf = accu(h);
    double sigma2_raw = accu(square(resid)) / std::max<double>(1, (x.n_elem - 1));
    double sigma2_edf = accu(square(resid)) / std::max<double>(1, (x.n_elem - std::max(1.0, edf)));
    vec w2 = weight_sq(x, grid, span, k_max);

    mat boot;
    if (keep_boot || adjust_mode == 1 || ci_type > 0)
        boot.set_size(grid.n_elem, B);
    vec below(grid.n_elem, fill::zeros);

    double alpha = 1 - level;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        std::mt19937_64 rng(0xB5297A4DUL
#ifdef _OPENMP
                            + (unsigned)omp_get_thread_num()
#endif
        );
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
        for (int b = 0; b < B; ++b)
        {
            // stratified residual resampling
            vec yb = fit_obs;
            for (uword k = 0; k < K; ++k)
            {
                const uvec &id = strata[k];
                if (id.n_elem == 0)
                    continue;
                for (uword t = 0; t < id.n_elem; ++t)
                {
                    uword pick = id[(uword)(rng() % id.n_elem)];
                    yb[id[t]] += resid[pick];
                }
            }
            vec fit_b = local_loess(x, yb, grid, span, deg, k_max);
            if (keep_boot || adjust_mode == 1 || ci_type > 0)
            {
                // Different b writes different columns, no data race
                boot.col(b) = fit_b;
            }
        }
    }

    // Quantiles helper
    auto qfun = [](vec v, double p) -> double
    {
        v = v(find_finite(v));
        if (v.n_elem == 0)
            return NA_REAL;
        v = sort(v);
        double h = p * (v.n_elem + 1.0);
        if (h <= 1.0)
            return v[0];
        if (h >= (double)v.n_elem)
            return v[v.n_elem - 1];
        double hf = std::floor(h);
        double dh = h - hf;
        uword i0 = (uword)hf - 1;
        uword i1 = std::min<uword>(v.n_elem - 1, i0 + 1);
        return v[i0] + dh * (v[i1] - v[i0]);
    };

    vec lo(grid.n_elem, fill::value(datum::nan)),
        hi(grid.n_elem, fill::value(datum::nan)),
        lo_basic(grid.n_elem, fill::value(datum::nan)),
        hi_basic(grid.n_elem, fill::value(datum::nan)),
        lo_bc(grid.n_elem, fill::value(datum::nan)),
        hi_bc(grid.n_elem, fill::value(datum::nan)),
        z0(grid.n_elem, fill::value(datum::nan));

    if (boot.n_elem)
    {
        for (uword j = 0; j < grid.n_elem; ++j)
        {
            vec bj = boot.row(j).t();
            double fq = qfun(bj, alpha / 2);
            double hq = qfun(bj, 1 - alpha / 2);
            lo[j] = fq;
            hi[j] = hq;
            // basic
            lo_basic[j] = 2 * fit_grid[j] - hq;
            hi_basic[j] = 2 * fit_grid[j] - fq;
            // BC (simple z0)
            uword cnt = accu(bj < fit_grid[j]);
            double prop = (cnt + 0.5) / (bj.n_elem + 1.0);
            double z = R::qnorm5(prop, 0, 1, 1, 0);
            z0[j] = z;
            // Simple BCa (without acceleration)
            double zalpha1 = R::qnorm5(alpha / 2, 0, 1, 1, 0);
            double zalpha2 = R::qnorm5(1 - alpha / 2, 0, 1, 1, 0);
            double qa = R::pnorm5(z0[j] + (z0[j] + zalpha1), 0, 1, 1, 0);
            double qb = R::pnorm5(z0[j] + (z0[j] + zalpha2), 0, 1, 1, 0);
            lo_bc[j] = qfun(bj, qa);
            hi_bc[j] = qfun(bj, qb);
        }
    }

    // Adjust mode (analytic scaling placeholder – here no extra scaling, kept for interface)
    if (adjust_mode == 1 && boot.n_elem)
    {
        // Users can still use lo/hi on R side; interface preserved without algorithm changes
    }

    // Residual MAD (robust scale)
    vec rc = resid - median(resid);
    double resid_mad = 1.4826 * median(abs(rc));

    List out = List::create(
        _["fit"] = fit_grid,
        _["lo"] = lo,
        _["hi"] = hi,
        _["lo_basic"] = lo_basic,
        _["hi_basic"] = hi_basic,
        _["lo_bc"] = lo_bc,
        _["hi_bc"] = hi_bc,
        _["z0"] = z0,
        _["below"] = below,
        _["B"] = B,
        _["edf"] = edf,
        _["sigma2_raw"] = sigma2_raw,
        _["sigma2_edf"] = sigma2_edf,
        _["w2"] = w2,
        _["alpha"] = alpha,
        _["level"] = level);
    if (keep_boot)
        out["boot"] = boot;
    out["x_obs"] = x;
    out["fit_obs"] = fit_obs;
    out["resid"] = resid;
    out["resid_mad"] = resid_mad;
    return out;
}
