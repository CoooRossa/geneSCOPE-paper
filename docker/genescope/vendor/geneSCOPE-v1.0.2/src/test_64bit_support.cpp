#include <Rcpp.h>

//' @title 检测当前编译/运行环境 64 位索引支持
//' @description 返回是否具备典型 64 位指针与 size_t (>=8 字节)，
//'   供在 R 端决定是否安全处理超大矩阵/稀疏结构（>2^31 行/列时需要）。
//' @return TRUE/FALSE
//' @export
// [[Rcpp::export]]
bool test_64bit_support()
{
    return (sizeof(void *) >= 8) && (sizeof(size_t) >= 8);
}
