#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
NumericVector objfcpp(NumericVector eps2, NumericVector sigma2, SEXP omega_,SEXP alpha_,SEXP beta_,SEXP r_,SEXP n_, SEXP tol_) 
{
	double omega = as<double>(omega_);
	double alpha = as<double>(alpha_);
	double beta = as<double>(beta_);
	double tol = as<double>(tol_);
	int r = as<int>(r_);
	int n = as<int>(n_);
	NumericVector vec(n);
	vec[0] = sigma2[0];
	for (int i = 1; i <= n-1; i++)
	{
		if (i > r-1)
		{
			vec[i] = max(tol, omega + alpha*eps2[i - 1] + beta*vec[i - 1]);
		}
		else
		{
			vec[i] = omega + alpha*eps2[i - 1] + beta*vec[i - 1];
		}
	}
	return vec;
}


// [[Rcpp::export]]

NumericVector objfcpp2(NumericVector eps2, NumericVector sigma2, SEXP omega_, SEXP alpha_, SEXP beta_, SEXP n_)
{
	double omega = as<double>(omega_);
	double alpha = as<double>(alpha_);
	double beta = as<double>(beta_);
	int n = as<int>(n_);
	NumericVector vec(n);
	vec[0] = sigma2[0];
	for (int i = 1; i <= n - 1; i++)
	{
		vec[i] = omega + alpha*eps2[i - 1] + beta*vec[i - 1];
	}
	return vec;
}
