// 1.SpatialWeight.cpp (2025-06-27)
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppEigen)]]
// Requires OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
using namespace Rcpp;

namespace {

inline bool enable_parallel_grid_nb() {
    const char* flag = std::getenv("GENESCOPE_ENABLE_OMP_GRID_NB");
    if (flag == nullptr) return false;
    if (flag[0] == '1') return true;
    if (flag[0] == 'T' || flag[0] == 't') return true;
    if (flag[0] == 'Y' || flag[0] == 'y') return true;
    if (flag[0] == 'F' || flag[0] == 'f') return true;
    if (flag[0] == 'U' || flag[0] == 'u') return true;
    return false;
}

struct GridGuardContext {
    const char* caller;
    int nrow;
    int ncol;
    long long n_cells;
};

inline GridGuardContext make_grid_context(const int nrow,
                                          const int ncol,
                                          const char* caller) {
    if (nrow == NA_INTEGER || ncol == NA_INTEGER) {
        Rcpp::stop("[%s] Grid dimensions must be finite integers.", caller);
    }
    if (nrow <= 0 || ncol <= 0) {
        Rcpp::stop("[%s] Grid dimensions must be positive (got %d x %d).",
                   caller, nrow, ncol);
    }
    const long long row_ll = static_cast<long long>(nrow);
    const long long col_ll = static_cast<long long>(ncol);
    const long long total = row_ll * col_ll;
    if (total > std::numeric_limits<int>::max()) {
        Rcpp::stop("[%s] Grid has %lld cells which exceeds INT_MAX (%d). Dimensions: %d x %d.",
                   caller, total, std::numeric_limits<int>::max(), nrow, ncol);
    }
    return GridGuardContext{caller, nrow, ncol, total};
}

inline long long safe_positive_mul(const GridGuardContext& ctx,
                                   const long long a,
                                   const long long b,
                                   const char* label) {
    if (a < 0 || b < 0) {
        Rcpp::stop("[%s] %s expects non-negative operands (grid %d x %d).",
                   ctx.caller, label, ctx.nrow, ctx.ncol);
    }
    if (a == 0 || b == 0) return 0;
    if (a > std::numeric_limits<long long>::max() / b) {
        Rcpp::stop("[%s] %s overflow for grid %d x %d (cells=%lld). Reduce grid resolution or adjust units.",
                   ctx.caller, label, ctx.nrow, ctx.ncol, ctx.n_cells);
    }
    return a * b;
}

inline long long safe_positive_add(const GridGuardContext& ctx,
                                   const long long a,
                                   const long long b,
                                   const char* label) {
    if (a < 0 || b < 0) {
        Rcpp::stop("[%s] %s expects non-negative operands (grid %d x %d).",
                   ctx.caller, label, ctx.nrow, ctx.ncol);
    }
    if (b > std::numeric_limits<long long>::max() - a) {
        Rcpp::stop("[%s] %s overflow for grid %d x %d (cells=%lld).",
                   ctx.caller, label, ctx.nrow, ctx.ncol, ctx.n_cells);
    }
    return a + b;
}

inline void validate_rect_edge_bound(const GridGuardContext& ctx,
                                     const bool queen) {
    const long long nrow_ll = static_cast<long long>(ctx.nrow);
    const long long ncol_ll = static_cast<long long>(ctx.ncol);
    const long long horizontal = safe_positive_mul(
        ctx, nrow_ll, std::max<long long>(ncol_ll - 1, 0), "horizontal edge count");
    const long long vertical = safe_positive_mul(
        ctx, ncol_ll, std::max<long long>(nrow_ll - 1, 0), "vertical edge count");
    long long undirected = safe_positive_add(ctx, horizontal, vertical, "rook adjacency total");
    if (queen) {
        const long long diag_pairs = safe_positive_mul(
            ctx, std::max<long long>(nrow_ll - 1, 0), std::max<long long>(ncol_ll - 1, 0), "diagonal adjacency count");
        const long long diag_edges = safe_positive_mul(ctx, diag_pairs, 2, "diagonal edge multiplier");
        undirected = safe_positive_add(ctx, undirected, diag_edges, "queen adjacency total");
    }
    safe_positive_mul(ctx, undirected, 2, "directed adjacency total");
}

inline void validate_hex_edge_bound(const GridGuardContext& ctx) {
    safe_positive_mul(ctx, ctx.n_cells, 6, "hex adjacency capacity");
}

inline List build_rectangular_nb(const int nrow,
                                 const int ncol,
                                 const bool queen,
                                 const bool allow_parallel) {
    const GridGuardContext ctx = make_grid_context(nrow, ncol, "grid_nb_impl");
    validate_rect_edge_bound(ctx, queen);
    const int n = static_cast<int>(ctx.n_cells);
    std::vector<std::vector<int>> nb(n);

    const int drow[8] = {0, 0, -1, 1, -1, -1, 1, 1};
    const int dcol[8] = {-1, 1, 0, 0, -1, 1, -1, 1};
    const int nOff = queen ? 8 : 4;

    const bool run_parallel =
#ifdef _OPENMP
        allow_parallel && enable_parallel_grid_nb();
#else
        false;
#endif

#ifdef _OPENMP
    if (run_parallel) {
#pragma omp parallel for schedule(static)
        for (int idx = 0; idx < n; ++idx) {
            int r = idx / ncol, c = idx % ncol;
            std::vector<int> local_nb;
            local_nb.reserve(nOff);
            for (int k = 0; k < nOff; ++k) {
                int rr = r + drow[k], cc = c + dcol[k];
                if (rr >= 0 && rr < nrow && cc >= 0 && cc < ncol) {
                    local_nb.push_back(rr * ncol + cc + 1);
                }
            }
            nb[idx].swap(local_nb);
        }
    } else
#endif
    {
        for (int idx = 0; idx < n; ++idx) {
            int r = idx / ncol, c = idx % ncol;
            std::vector<int> local_nb;
            local_nb.reserve(nOff);
            for (int k = 0; k < nOff; ++k) {
                int rr = r + drow[k], cc = c + dcol[k];
                if (rr >= 0 && rr < nrow && cc >= 0 && cc < ncol) {
                    local_nb.push_back(rr * ncol + cc + 1);
                }
            }
            nb[idx].swap(local_nb);
        }
    }

    List out(n);
    for (int i = 0; i < n; ++i)
        out[i] = IntegerVector(nb[i].begin(), nb[i].end());
    out.attr("class") = "nb";
    out.attr("region.id") = Rcpp::seq(1, n);
    out.attr("queen") = queen;
    return out;
}

inline List build_hex_r_nb(const int nrow,
                           const int ncol,
                           const bool oddr,
                           const bool allow_parallel) {
    const GridGuardContext ctx = make_grid_context(nrow, ncol, "grid_nb_hex_impl");
    validate_hex_edge_bound(ctx);
    const int n = static_cast<int>(ctx.n_cells);
    std::vector<std::vector<int>> nb(n);

    const bool run_parallel =
#ifdef _OPENMP
        allow_parallel && enable_parallel_grid_nb();
#else
        false;
#endif

#ifdef _OPENMP
    if (run_parallel) {
#pragma omp parallel for schedule(static)
        for (int idx = 0; idx < n; ++idx) {
            int r = idx / ncol, c = idx % ncol;
            std::vector<int> local_nb;
            local_nb.reserve(6);
            const bool is_odd = (r % 2) == 1;
            int rr = r, cc = c - 1;
            if (cc >= 0) local_nb.push_back(rr * ncol + cc + 1);
            cc = c + 1;
            if (cc < ncol) local_nb.push_back(rr * ncol + cc + 1);
            if (oddr) {
                if (r - 1 >= 0) {
                    cc = c + (is_odd ? 1 : 0);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                    cc = c + (is_odd ? 0 : -1);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                }
                if (r + 1 < nrow) {
                    cc = c + (is_odd ? 1 : 0);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                    cc = c + (is_odd ? 0 : -1);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                }
            } else {
                if (r - 1 >= 0) {
                    cc = c + (is_odd ? 0 : 1);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                    cc = c + (is_odd ? -1 : 0);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                }
                if (r + 1 < nrow) {
                    cc = c + (is_odd ? 0 : 1);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                    cc = c + (is_odd ? -1 : 0);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                }
            }
            nb[idx].swap(local_nb);
        }
    } else
#endif
    {
        for (int idx = 0; idx < n; ++idx) {
            int r = idx / ncol, c = idx % ncol;
            std::vector<int> local_nb;
            local_nb.reserve(6);
            const bool is_odd = (r % 2) == 1;
            int rr = r, cc = c - 1;
            if (cc >= 0) local_nb.push_back(rr * ncol + cc + 1);
            cc = c + 1;
            if (cc < ncol) local_nb.push_back(rr * ncol + cc + 1);
            if (oddr) {
                if (r - 1 >= 0) {
                    cc = c + (is_odd ? 1 : 0);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                    cc = c + (is_odd ? 0 : -1);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                }
                if (r + 1 < nrow) {
                    cc = c + (is_odd ? 1 : 0);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                    cc = c + (is_odd ? 0 : -1);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                }
            } else {
                if (r - 1 >= 0) {
                    cc = c + (is_odd ? 0 : 1);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                    cc = c + (is_odd ? -1 : 0);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r - 1) * ncol + cc + 1);
                }
                if (r + 1 < nrow) {
                    cc = c + (is_odd ? 0 : 1);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                    cc = c + (is_odd ? -1 : 0);
                    if (cc >= 0 && cc < ncol) local_nb.push_back((r + 1) * ncol + cc + 1);
                }
            }
            nb[idx].swap(local_nb);
        }
    }

    List out(n);
    for (int i = 0; i < n; ++i)
        out[i] = IntegerVector(nb[i].begin(), nb[i].end());
    out.attr("class") = "nb";
    out.attr("region.id") = Rcpp::seq(1, n);
    out.attr("queen") = false;
    out.attr("topology") = oddr ? "hex-oddr" : "hex-evenr";
    return out;
}

inline List build_hex_q_nb(const int nrow,
                           const int ncol,
                           const bool oddq,
                           const bool allow_parallel) {
    List nb_t = build_hex_r_nb(ncol, nrow, oddq, allow_parallel);
    nb_t.attr("queen") = false;
    nb_t.attr("topology") = oddq ? "hex-oddq" : "hex-evenq";
    return nb_t;
}

}

/*------------------------------------------------------------------*/
/* 1. Parallel generation of Queen/Rook adjacency list            */
/*------------------------------------------------------------------*/
//' @title Build Queen or Rook Neighbourhood List in Parallel
//'
//' @description
//'   Creates a spatial neighbourhood list for a rectangular grid using either
//'   Queen (8-neighbour) or Rook (4-neighbour) connectivity. OpenMP
//'   parallelisation scales efficiently to millions of grid cells.
//'
//' @param nrow Integer. Number of rows in the grid.
//' @param ncol Integer. Number of columns in the grid.
//' @param queen Logical. TRUE (default) for Queen moves, FALSE for Rook.
//'
//' @return A list of integer vectors of class \code{nb}. Indices are 1-based
//'   and attributes \code{region.id} and \code{queen} are set.
//'
//' @details
//'   Offsets for rook and (optionally) diagonal neighbours are pre-computed and
//'   applied in an OpenMP \code{for} loop. Each thread accumulates a private
//'   neighbour vector which is swapped into the main container to avoid
//'   locking.
//'
//' @examples
//' \dontrun{
//' nb <- grid_nb_omp(1024, 1024, queen = TRUE)
//' head(nb[[1]])
//' }
// [[Rcpp::export]]
List grid_nb_omp(const int nrow,
                 const int ncol,
                 const bool queen = true)
{
    return build_rectangular_nb(nrow, ncol, queen, true);
}

// [[Rcpp::export]]
List grid_nb_serial(const int nrow,
                    const int ncol,
                    const bool queen = true)
{
    return build_rectangular_nb(nrow, ncol, queen, false);
}

/*------------------------------------------------------------------*/
/* 1b. Parallel generation of Hex adjacency (odd-r / even-r)       */
/*------------------------------------------------------------------*/
//' @title Build Hexagonal (6-neighbour) Neighbourhood List in Parallel
//'
//' @description
//'   Creates a spatial neighbourhood list for a hex tiling laid out on a
//'   rectangular index with row-wise offsets (odd-r or even-r). See
//'   https://www.redblobgames.com/grids/hex-grids/ for conventions.
//'
//' @param nrow Integer. Number of rows in the logical lattice.
//' @param ncol Integer. Number of columns in the logical lattice.
//' @param oddr Logical. If TRUE, use odd-r layout (odd rows shifted right);
//'   otherwise use even-r layout (even rows shifted right).
//'
//' @return A list of integer vectors of class \code{nb}. Indices are 1-based
//'   and attributes \code{region.id} and \code{topology} ("hex-oddr" or
//'   "hex-evenr") are set. Attribute \code{queen} is set to FALSE for
//'   compatibility with existing code expecting a boolean.
//'
//' @examples
//' \dontrun{
//' nb_hex <- grid_nb_hex_omp(256, 512, oddr = TRUE)
//' }
// [[Rcpp::export]]
List grid_nb_hex_omp(const int nrow,
                     const int ncol,
                     const bool oddr = true)
{
    return build_hex_r_nb(nrow, ncol, oddr, true);
}

// [[Rcpp::export]]
List grid_nb_hex(const int nrow,
                 const int ncol,
                 const bool oddr = true)
{
    return build_hex_r_nb(nrow, ncol, oddr, false);
}

/*------------------------------------------------------------------*/
/* 1c. Parallel generation of Hex adjacency (odd-q / even-q)       */
/*------------------------------------------------------------------*/
//' @title Build Hexagonal (6-neighbour) Neighbour List with column-offset (odd-q/even-q)
//'
//' @description
//'   Convenience wrapper that reuses the row-offset builder on a transposed
//'   lattice to implement pointy-top hex coordinates with column offsets
//'   (odd-q/even-q). Region indices are 1..(nrow*ncol) in column-major order
//'   consistent with mapping id = (gx-1)*nrow + gy.
//'
//' @param nrow Integer. Number of rows in the logical lattice.
//' @param ncol Integer. Number of columns in the logical lattice.
//' @param oddq Logical. If TRUE, odd columns are shifted; otherwise even.
//'
//' @return An \code{nb} list with attributes \code{region.id}, \code{queen}=FALSE,
//'   and \code{topology} set to \code{"hex-oddq"} or \code{"hex-evenq"}.
// [[Rcpp::export]]
List grid_nb_hexq_omp(const int nrow,
                      const int ncol,
                      const bool oddq = true)
{
    return build_hex_q_nb(nrow, ncol, oddq, true);
}

// [[Rcpp::export]]
List grid_nb_hexq(const int nrow,
                  const int ncol,
                  const bool oddq = true)
{
    return build_hex_q_nb(nrow, ncol, oddq, false);
}

/*------------------------------------------------------------------*/
/* 2. Build binary sparse matrix (dgCMatrix) in parallel          */
/*------------------------------------------------------------------*/
//' @title Convert Neighbour List to Binary Sparse Matrix (dgCMatrix)
//'
//' @description
//'   Transforms an \code{nb} neighbourhood list into a column-compressed
//'   binary sparse matrix of class \code{dgCMatrix}. The routine streams the
//'   triplet representation in parallel with OpenMP and finally constructs the
//'   sparse matrix with Eigen.
//'
//' @param nb An object of class \code{nb}, typically returned by
//'   \code{grid_nb_omp()} or \pkg{spdep} functions.
//'
//' @return A \code{dgCMatrix} of dimension \eqn{n \times n} where \eqn{n} is
//'   the number of regions. Entries are \code{1} where two regions are
//'   adjacent and \code{0} elsewhere.
//'
//' @details
//'   The neighbour indices are 1-based; they are converted to 0-based before
//'   creation of the triplets for Eigen. Each thread accumulates its own
//'   vector of triplets to avoid synchronisation. After the parallel region
//'   all triplets are concatenated and the resulting sparse matrix is
//'   compressed in column-major order.
//'
//' @examples
//' \dontrun{
//' nb  <- grid_nb_omp(512, 512, queen = FALSE)
//' W   <- listw_B_omp(nb)
//' Matrix::nnzero(W) / length(nb)  # average degree
//' }
// [[Rcpp::export]]
SEXP listw_B_omp(const List nb)
{
    const int n = nb.size();

    /* 2.1 Single-thread deep copy into C++ containers — fully isolate from R memory */
    std::vector<std::vector<int>> nb_cpp(n);
    for (int i = 0; i < n; ++i)
    {
        IntegerVector v = nb[i];
        nb_cpp[i].assign(v.begin(), v.end());
    }

    /* 2.2 Thread-private triplet buffers */
#ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#else
    const int nThreads = 1;
#endif
    std::vector<std::vector<Eigen::Triplet<double>>> trip_tls(nThreads);

#pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        int tid =
#ifdef _OPENMP
            omp_get_thread_num();
#else
            0;
#endif
        auto &vec = trip_tls[tid];
        vec.reserve(vec.size() + nb_cpp[i].size());
        for (int j : nb_cpp[i])
        {
            if (j == 0)
                continue;
            vec.emplace_back(i, j - 1, 1.0); // 0-based
        }
    }

    /* 2.3 Merge triplets */
    size_t total = 0;
    for (auto &v : trip_tls)
        total += v.size();
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(total);
    for (auto &v : trip_tls)
    {
        std::move(v.begin(), v.end(), std::back_inserter(trip));
        std::vector<Eigen::Triplet<double>>().swap(v); // release thread-local buffer
    }

    /* 2.4 Construct compressed sparse matrix */
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> W(n, n);
    W.setFromTriplets(trip.begin(), trip.end());
    W.makeCompressed();

    return Rcpp::wrap(W); // dgCMatrix
}

/*------------------------------------------------------------------*/
/* 3. Kernel-smoothed sparse weights for active cells              */
/*------------------------------------------------------------------*/

namespace {

struct RectOffset {
    int dr;
    int dc;
    double dist;
};

inline std::vector<RectOffset> make_rect_offsets(const int radius, const bool queen) {
    std::vector<RectOffset> offsets;
    if (radius < 1) return offsets;
    offsets.reserve((2 * radius + 1) * (2 * radius + 1) - 1);
    for (int dr = -radius; dr <= radius; ++dr) {
        for (int dc = -radius; dc <= radius; ++dc) {
            if (dr == 0 && dc == 0) continue;
            if (!queen && (std::abs(dr) + std::abs(dc) > radius)) continue;
            double d = std::sqrt(static_cast<double>(dr * dr + dc * dc));
            offsets.push_back({dr, dc, d});
        }
    }
    return offsets;
}

struct CubeOffset {
    int dx;
    int dy;
    int dz;
    int dist;
};

inline std::vector<CubeOffset> make_cube_offsets(const int radius) {
    std::vector<CubeOffset> out;
    if (radius < 1) return out;
    out.reserve(3 * radius * (radius + 1)); // exclude self
    for (int dx = -radius; dx <= radius; ++dx) {
        int dy_min = std::max(-radius, -dx - radius);
        int dy_max = std::min(radius, -dx + radius);
        for (int dy = dy_min; dy <= dy_max; ++dy) {
            int dz = -dx - dy;
            if (dx == 0 && dy == 0 && dz == 0) continue;
            int dist = std::max({std::abs(dx), std::abs(dy), std::abs(dz)});
            out.push_back({dx, dy, dz, dist});
        }
    }
    return out;
}

inline double raw_weight(const int dist, const bool flat, const double sigma) {
    if (flat) return 1.0;
    const double s2 = sigma * sigma;
    return std::exp(-(static_cast<double>(dist * dist)) / (2.0 * s2));
}

// RedBlobGames conversions: odd-r / even-r (row-offset)
inline std::pair<int, int> offset_to_axial_oddr(const int row, const int col) {
    int q = col - ((row - (row & 1)) / 2);
    int r = row;
    return {q, r};
}

inline std::pair<int, int> offset_to_axial_evenr(const int row, const int col) {
    int q = col - ((row + (row & 1)) / 2);
    int r = row;
    return {q, r};
}

inline std::pair<int, int> axial_to_offset_oddr(const int q, const int r) {
    int col = q + ((r - (r & 1)) / 2);
    int row = r;
    return {row, col};
}

inline std::pair<int, int> axial_to_offset_evenr(const int q, const int r) {
    int col = q + ((r + (r & 1)) / 2);
    int row = r;
    return {row, col};
}

// RedBlobGames conversions: odd-q / even-q (col-offset)
inline std::pair<int, int> offset_to_axial_oddq(const int row, const int col) {
    int q = col;
    int r = row - ((col - (col & 1)) / 2);
    return {q, r};
}

inline std::pair<int, int> offset_to_axial_evenq(const int row, const int col) {
    int q = col;
    int r = row - ((col + (col & 1)) / 2);
    return {q, r};
}

inline std::pair<int, int> axial_to_offset_oddq(const int q, const int r) {
    int col = q;
    int row = r + ((q - (q & 1)) / 2);
    return {row, col};
}

inline std::pair<int, int> axial_to_offset_evenq(const int q, const int r) {
    int col = q;
    int row = r + ((q + (q & 1)) / 2);
    return {row, col};
}

} // namespace

// [[Rcpp::export]]
SEXP grid_weights_kernel_rect_omp(
    const int nrow,
    const int ncol,
    const IntegerVector gx,   // 1-based active x
    const IntegerVector gy,   // 1-based active y
    const bool queen = true,
    const int radius = 2,
    const std::string kernel = "gaussian", // "gaussian" | "flat"
    const double sigma = 1.0
) {
    const int n_active = gx.size();
    if (gy.size() != n_active) {
        stop("gx and gy must have the same length.");
    }
    if (nrow <= 0 || ncol <= 0) {
        stop("nrow and ncol must be positive.");
    }
    const int rad = radius;
    if (rad < 1) {
        stop("radius must be >= 1.");
    }
    const bool flat = (kernel == "flat");
    if (!flat && kernel != "gaussian") {
        stop("kernel must be 'gaussian' or 'flat'.");
    }
    const double sigma_use = (sigma > 0.0) ? sigma : 1.0;

    const int n_full = nrow * ncol;
    std::vector<int> idx_map(n_full, -1);
    std::vector<int> active_row(n_active);
    std::vector<int> active_col(n_active);

    for (int i = 0; i < n_active; ++i) {
        int col0 = gx[i] - 1;
        int row0 = gy[i] - 1;
        if (col0 < 0 || col0 >= ncol || row0 < 0 || row0 >= nrow) {
            stop("gx/gy out of bounds for the provided nrow/ncol.");
        }
        active_row[i] = row0;
        active_col[i] = col0;
        int full_id0 = row0 * ncol + col0; // row-major
        idx_map[full_id0] = i;
    }

    const std::vector<RectOffset> offsets = make_rect_offsets(rad, queen);

#ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#else
    const int nThreads = 1;
#endif
    std::vector<std::vector<Eigen::Triplet<double>>> trip_tls(nThreads);

#pragma omp parallel for schedule(static)
    for (int i_act = 0; i_act < n_active; ++i_act) {
        int tid =
#ifdef _OPENMP
            omp_get_thread_num();
#else
            0;
#endif
        auto &trip = trip_tls[tid];

        const int r0 = active_row[i_act];
        const int c0 = active_col[i_act];

        std::vector<int> neigh_idx;
        std::vector<double> neigh_raw;
        neigh_idx.reserve(offsets.size());
        neigh_raw.reserve(offsets.size());
        double wsum = 0.0;

        for (const auto &off : offsets) {
            int rr = r0 + off.dr;
            int cc = c0 + off.dc;
            if (rr < 0 || rr >= nrow || cc < 0 || cc >= ncol) continue;
            int full_j0 = rr * ncol + cc;
            int j_act = idx_map[full_j0];
            if (j_act < 0) continue; // inactive
            double w_raw = flat ? 1.0 : std::exp(-(off.dist * off.dist) / (2.0 * sigma_use * sigma_use));
            neigh_idx.push_back(j_act);
            neigh_raw.push_back(w_raw);
            wsum += w_raw;
        }

        if (wsum > 0.0) {
            trip.reserve(trip.size() + neigh_idx.size());
            for (size_t k = 0; k < neigh_idx.size(); ++k) {
                trip.emplace_back(i_act, neigh_idx[k], neigh_raw[k] / wsum);
            }
        }
    }

    size_t total = 0;
    for (auto &v : trip_tls) total += v.size();
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(total);
    for (auto &v : trip_tls) {
        std::move(v.begin(), v.end(), std::back_inserter(trip));
        std::vector<Eigen::Triplet<double>>().swap(v);
    }

    Eigen::SparseMatrix<double, Eigen::ColMajor, int> W(n_active, n_active);
    W.setFromTriplets(trip.begin(), trip.end());
    W.makeCompressed();
    return Rcpp::wrap(W);
}

// [[Rcpp::export]]
SEXP grid_weights_kernel_hexr_omp(
    const int nrow,
    const int ncol,
    const IntegerVector gx,   // 1-based active x
    const IntegerVector gy,   // 1-based active y
    const bool oddr = true,
    const int radius = 2,
    const std::string kernel = "gaussian",
    const double sigma = 1.0
) {
    const int n_active = gx.size();
    if (gy.size() != n_active) stop("gx and gy must have the same length.");
    if (nrow <= 0 || ncol <= 0) stop("nrow and ncol must be positive.");
    const int rad = radius;
    if (rad < 1) stop("radius must be >= 1.");
    const bool flat = (kernel == "flat");
    if (!flat && kernel != "gaussian") stop("kernel must be 'gaussian' or 'flat'.");
    const double sigma_use = (sigma > 0.0) ? sigma : 1.0;

    const int n_full = nrow * ncol;
    std::vector<int> idx_map(n_full, -1);
    std::vector<int> active_row(n_active);
    std::vector<int> active_col(n_active);

    for (int i = 0; i < n_active; ++i) {
        int col0 = gx[i] - 1;
        int row0 = gy[i] - 1;
        if (col0 < 0 || col0 >= ncol || row0 < 0 || row0 >= nrow) {
            stop("gx/gy out of bounds for the provided nrow/ncol.");
        }
        active_row[i] = row0;
        active_col[i] = col0;
        int full_id0 = row0 * ncol + col0; // row-major (matches grid_nb_hex_omp)
        idx_map[full_id0] = i;
    }

    const std::vector<CubeOffset> cube_offs = make_cube_offsets(rad);

#ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#else
    const int nThreads = 1;
#endif
    std::vector<std::vector<Eigen::Triplet<double>>> trip_tls(nThreads);

#pragma omp parallel for schedule(static)
    for (int i_act = 0; i_act < n_active; ++i_act) {
        int tid =
#ifdef _OPENMP
            omp_get_thread_num();
#else
            0;
#endif
        auto &trip = trip_tls[tid];

        const int r0 = active_row[i_act];
        const int c0 = active_col[i_act];
        std::pair<int, int> axial0 = oddr ? offset_to_axial_oddr(r0, c0) : offset_to_axial_evenr(r0, c0);
        const int q0 = axial0.first;
        const int ar0 = axial0.second;

        const int x0 = q0;
        const int z0 = ar0;
        const int y0 = -x0 - z0;

        std::vector<int> neigh_idx;
        std::vector<double> neigh_raw;
        neigh_idx.reserve(cube_offs.size());
        neigh_raw.reserve(cube_offs.size());
        double wsum = 0.0;

        for (const auto &off : cube_offs) {
            int x1 = x0 + off.dx;
            int y1 = y0 + off.dy;
            int z1 = z0 + off.dz;
            int q1 = x1;
            int r1 = z1;
            std::pair<int, int> offset1 = oddr ? axial_to_offset_oddr(q1, r1) : axial_to_offset_evenr(q1, r1);
            int rr = offset1.first;
            int cc = offset1.second;
            if (rr < 0 || rr >= nrow || cc < 0 || cc >= ncol) continue;
            int full_j0 = rr * ncol + cc;
            int j_act = idx_map[full_j0];
            if (j_act < 0) continue;
            double w_raw = raw_weight(off.dist, flat, sigma_use);
            neigh_idx.push_back(j_act);
            neigh_raw.push_back(w_raw);
            wsum += w_raw;
        }

        if (wsum > 0.0) {
            trip.reserve(trip.size() + neigh_idx.size());
            for (size_t k = 0; k < neigh_idx.size(); ++k) {
                trip.emplace_back(i_act, neigh_idx[k], neigh_raw[k] / wsum);
            }
        }
    }

    size_t total = 0;
    for (auto &v : trip_tls) total += v.size();
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(total);
    for (auto &v : trip_tls) {
        std::move(v.begin(), v.end(), std::back_inserter(trip));
        std::vector<Eigen::Triplet<double>>().swap(v);
    }

    Eigen::SparseMatrix<double, Eigen::ColMajor, int> W(n_active, n_active);
    W.setFromTriplets(trip.begin(), trip.end());
    W.makeCompressed();
    return Rcpp::wrap(W);
}

// [[Rcpp::export]]
SEXP grid_weights_kernel_hexq_omp(
    const int nrow,
    const int ncol,
    const IntegerVector gx,   // 1-based active x
    const IntegerVector gy,   // 1-based active y
    const bool oddq = true,
    const int radius = 2,
    const std::string kernel = "gaussian",
    const double sigma = 1.0
) {
    const int n_active = gx.size();
    if (gy.size() != n_active) stop("gx and gy must have the same length.");
    if (nrow <= 0 || ncol <= 0) stop("nrow and ncol must be positive.");
    const int rad = radius;
    if (rad < 1) stop("radius must be >= 1.");
    const bool flat = (kernel == "flat");
    if (!flat && kernel != "gaussian") stop("kernel must be 'gaussian' or 'flat'.");
    const double sigma_use = (sigma > 0.0) ? sigma : 1.0;

    const int n_full = nrow * ncol;
    std::vector<int> idx_map(n_full, -1);
    std::vector<int> active_row(n_active);
    std::vector<int> active_col(n_active);

    for (int i = 0; i < n_active; ++i) {
        int col0 = gx[i] - 1;
        int row0 = gy[i] - 1;
        if (col0 < 0 || col0 >= ncol || row0 < 0 || row0 >= nrow) {
            stop("gx/gy out of bounds for the provided nrow/ncol.");
        }
        active_row[i] = row0;
        active_col[i] = col0;
        int full_id0 = col0 * nrow + row0; // column-major (matches grid_nb_hexq_omp)
        idx_map[full_id0] = i;
    }

    const std::vector<CubeOffset> cube_offs = make_cube_offsets(rad);

#ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#else
    const int nThreads = 1;
#endif
    std::vector<std::vector<Eigen::Triplet<double>>> trip_tls(nThreads);

#pragma omp parallel for schedule(static)
    for (int i_act = 0; i_act < n_active; ++i_act) {
        int tid =
#ifdef _OPENMP
            omp_get_thread_num();
#else
            0;
#endif
        auto &trip = trip_tls[tid];

        const int r0 = active_row[i_act];
        const int c0 = active_col[i_act];
        std::pair<int, int> axial0 = oddq ? offset_to_axial_oddq(r0, c0) : offset_to_axial_evenq(r0, c0);
        const int q0 = axial0.first;
        const int ar0 = axial0.second;

        const int x0 = q0;
        const int z0 = ar0;
        const int y0 = -x0 - z0;

        std::vector<int> neigh_idx;
        std::vector<double> neigh_raw;
        neigh_idx.reserve(cube_offs.size());
        neigh_raw.reserve(cube_offs.size());
        double wsum = 0.0;

        for (const auto &off : cube_offs) {
            int x1 = x0 + off.dx;
            int y1 = y0 + off.dy;
            int z1 = z0 + off.dz;
            int q1 = x1;
            int r1 = z1;
            std::pair<int, int> offset1 = oddq ? axial_to_offset_oddq(q1, r1) : axial_to_offset_evenq(q1, r1);
            int rr = offset1.first;
            int cc = offset1.second;
            if (rr < 0 || rr >= nrow || cc < 0 || cc >= ncol) continue;
            int full_j0 = cc * nrow + rr; // column-major
            int j_act = idx_map[full_j0];
            if (j_act < 0) continue;
            double w_raw = raw_weight(off.dist, flat, sigma_use);
            neigh_idx.push_back(j_act);
            neigh_raw.push_back(w_raw);
            wsum += w_raw;
        }

        if (wsum > 0.0) {
            trip.reserve(trip.size() + neigh_idx.size());
            for (size_t k = 0; k < neigh_idx.size(); ++k) {
                trip.emplace_back(i_act, neigh_idx[k], neigh_raw[k] / wsum);
            }
        }
    }

    size_t total = 0;
    for (auto &v : trip_tls) total += v.size();
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(total);
    for (auto &v : trip_tls) {
        std::move(v.begin(), v.end(), std::back_inserter(trip));
        std::vector<Eigen::Triplet<double>>().swap(v);
    }

    Eigen::SparseMatrix<double, Eigen::ColMajor, int> W(n_active, n_active);
    W.setFromTriplets(trip.begin(), trip.end());
    W.makeCompressed();
    return Rcpp::wrap(W);
}

// Compatibility alias that calls the original function
// [[Rcpp::export]]
SEXP grid_nb(const int nrow,
             const int ncol,
             const bool queen = true)
{
    return grid_nb_serial(nrow, ncol, queen);
}

// [[Rcpp::export]]
SEXP nb2mat(SEXP nb)
{
    return listw_B_omp(nb);
}
