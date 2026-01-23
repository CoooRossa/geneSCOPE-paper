#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

Rcpp::List native_openmp_info() {
  Rcpp::List out;
#ifdef _OPENMP
  out["compiled_with_openmp"] = true;
  out["openmp_version"] = static_cast<int>(_OPENMP);
  out["omp_num_procs"] = omp_get_num_procs();
  out["omp_max_threads"] = omp_get_max_threads();
  out["omp_dynamic"] = static_cast<bool>(omp_get_dynamic());
  out["omp_nested"] = static_cast<bool>(omp_get_max_active_levels() > 1);
#else
  out["compiled_with_openmp"] = false;
  out["openmp_version"] = Rcpp::IntegerVector::create(NA_INTEGER);
  out["omp_num_procs"] = Rcpp::IntegerVector::create(NA_INTEGER);
  out["omp_max_threads"] = Rcpp::IntegerVector::create(NA_INTEGER);
  out["omp_dynamic"] = Rcpp::LogicalVector::create(NA_LOGICAL);
  out["omp_nested"] = Rcpp::LogicalVector::create(NA_LOGICAL);
#endif
  return out;
}

Rcpp::List native_openmp_set_threads(const int n_threads) {
  Rcpp::List out;
  out["requested_threads"] = n_threads;
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
  out["compiled_with_openmp"] = true;
  out["omp_max_threads"] = omp_get_max_threads();
  out["omp_num_procs"] = omp_get_num_procs();
#else
  out["compiled_with_openmp"] = false;
  out["omp_max_threads"] = Rcpp::IntegerVector::create(NA_INTEGER);
  out["omp_num_procs"] = Rcpp::IntegerVector::create(NA_INTEGER);
  out["note"] = "OpenMP not enabled at compile time (_OPENMP undefined).";
#endif
  return out;
}

RcppExport SEXP _geneSCOPERebuild_native_openmp_info(void) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  rcpp_result_gen = Rcpp::wrap(native_openmp_info());
  return rcpp_result_gen;
  END_RCPP
}

RcppExport SEXP _geneSCOPERebuild_native_openmp_set_threads(SEXP n_threadsSEXP) {
  BEGIN_RCPP
  Rcpp::RObject rcpp_result_gen;
  Rcpp::traits::input_parameter< const int >::type n_threads(n_threadsSEXP);
  rcpp_result_gen = Rcpp::wrap(native_openmp_set_threads(n_threads));
  return rcpp_result_gen;
  END_RCPP
}
