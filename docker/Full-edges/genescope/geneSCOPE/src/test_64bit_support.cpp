#include <Rcpp.h>

//' @title Detect 64-bit index support in the current build/runtime environment
//' @description Returns whether typical 64-bit pointers and size_t (>=8 bytes) are available,
//'   so the R side can decide whether it is safe to handle very large matrices/sparse structures
//'   (needed when >2^31 rows/cols).
//' @return TRUE/FALSE
//' @export
// [[Rcpp::export]]
bool test_64bit_support()
{
    return (sizeof(void *) >= 8) && (sizeof(size_t) >= 8);
}
