#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix cbind_cpp(NumericMatrix x, NumericMatrix y)
{
  unsigned long int x_nrow = x.nrow();
  unsigned long int x_ncol = x.ncol();
  unsigned long int y_nrow = y.nrow();
  unsigned long int y_ncol = y.ncol();
  if (x_nrow != y_nrow) {
    ::Rf_error("you must have x.nrow() == y.nrow()");
  }
  unsigned long int res_ncol = x_ncol + y_ncol;
  NumericMatrix res(x_nrow, res_ncol);

  for(int j = 0; j < res_ncol; j++) {
    if(j < x_ncol)
    {
      res(_ , j) = x(_ , j);
    } else {
      res(_, j) = y(_ , j - x_ncol);
    }
  }

  return res;
}
