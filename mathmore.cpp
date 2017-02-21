#include "mathmore.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>

long long binom(unsigned n, unsigned k)
{
    if (0 == k || n == k) 
        return 1;
    if (k > n) 
        return 0;
    if (k > (n - k)) 
        k = n - k;
    if (1 == k) 
        return n;
    long long b = 1;
    for (unsigned i = 1; i <= k; ++i)
    {
        if ( b > std::numeric_limits<long long>::max() / n-(k-i) )
        {
            std::cerr<<"binomial coefficients overflow!"<<std::endl;
            exit(1); /* Overflow */
        }
        b *= (n - (k - i));
        b /= i;
    }
    return b;
}


double mean( const std::vector<int> &v)
{
    return (std::accumulate(v.begin(), v.end(), 0.0) * 1.0 / v.size()); 
}

double rmse( const std::vector<int> &v, int truth )
{
    double sum = 0;
    for ( const auto &x : v )
        sum += 1.0*(x-truth)*(x-truth);
    return sqrt( sum / v.size() );
}

double stdev( const std::vector<int> &v )
{
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();

    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), std::bind2nd(std::minus<double>(), mean));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    return std::sqrt(sq_sum / (v.size()-1) );
}


double min_positive_normalized( const std::vector<double> &p ) // p may not be normalized.
{
    double sum = std::accumulate(p.begin(), p.end(), 0.0);
    double min = 1.0;
    for ( const auto &x : p )
        if ( (x>0) && (x/sum < min) )
            min = x/sum;
    return min;
}

double cnt_positive( const std::vector<double> &p )
{
    int cnt = 0;
    for ( const auto &x : p )
        if ( x>0 )
            cnt++;
    return cnt;
}













std::vector<double> ChebMore::expand_basic() const // Expansion of cos(L arccos x)
{
    std::vector<double> a(L+1);
	
    if (L % 2 == 0)
    {
        unsigned n = L / 2;
        for (unsigned k = 0; k <= n; k++)
            a[2*k] = pow(-1,n-k) * pow(2,2*k) * n / (n+k) * binom(n+k,2*k); // may be simpler
    }
    else
    {
        unsigned n = (L-1)/2;
        for (unsigned k = 0; k <= n; k++)
            a[2*k+1] = pow(-1,n-k) * pow(2,2*k) * (2*n+1) / (2*k+1) * binom(n+k,2*k); // may be simpler
    }
    return a;
}

std::vector<double> ChebMore::expand() const // Expansion of cos(L arccos (bx+c))
{
    std::vector<double> a0 = expand_basic();
    std::vector<double> a(L+1);
	
    // sum_i a[i] x^i = sum_i a0[i]*(bx+c)^i
    for (int i = 0; i <= L; i++)
        for (int j = 0; j <= i; j++)
            a[j] += a0[i] * pow(b,j) * pow(c,i-j) * binom(i,j);
    return a;
}

double ChebMore::evaluate(double x) const
{
    double A = b * x + c;
    if ( A > 1 )
    {
        return cosh(L*acosh(A));
    }
    else if ( A >= -1 )
    {
        return cos(L*acos(A));
    }
    else
    {
        return ( cosh(L*acosh(-A)) * ( (L%2==0)? 1:-1)) ;
    }
}
