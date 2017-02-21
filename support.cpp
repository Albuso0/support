#include "support.h"
#include "mathmore.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>




/* ----------------------POLYNOMIAL ESTIMATOR NO SPLITTING-------------------------------- */
double Support::estimate() const
{
    if ( pmin < Ratio / n ) // const
    {
        double result = 0.0;
        for ( const auto & pair : fin )
                result += getCoeff(pair.first) * pair.second;
        return result;
    }
    else
        return estimate_plug();
}

double Support::getCoeff( int N ) const
{
    if ( N == 0 ) return 0;
    if ( N > L ) return 1;

    double A = ( Ratio + n*pmin ) / ( Ratio - n*pmin );
    ChebMore cheb(L, 1, -A); // polynomial of cos L arccos(t-A)
    std::vector<double> a = cheb.expand(); // Expand: cos L arccos(x-A)=sum_i a[i]*x^i
    double amp = cheb.evaluate(0); // cos L arccos(-A)

    double result = -a[N]/amp;
    double s = 2 / ( Ratio - n * pmin );
    for ( int i = 1; i <= N; ++i )
        result *= i * s;

    return (1+result);
}
/* ----------------------END POLYNOMIAL ESTIMATOR NO SPLITTING-------------------------------- */





/* ----------------------OTHER ESTIMATORS-------------------------------- */
// Ref for TG, CL1, CL2: "Estimating the Number of Classes via Sample Coverage, Chao-Lee, 1992"
// Ref for JK: "Robust Estimation of Population Size When Capture Probabilities Vary Among Animals. K. P. Burnham and W. S. Overton, 1979"

double Support::estimate_plug() const
{
    double result = 0.0;
    for ( const auto & pair : fin )
        result += pair.second;
    return result;
}

double Support::coverage_TG() const
{
    double n1 = 0;
    for ( const auto & pair : fin )
        if ( pair.first == 1 )
        {
            n1 = pair.second;
            break;
        }
    return (1 - n1/n);
}

double Support::estimate_TG() const
{
    return ( estimate_plug() / coverage_TG() );
}


double Support::estimate_J1() const
{
    double f1 = 0;
    for ( const auto & pair : fin )
        if ( pair.first == 1 )
        {
            f1 = pair.second;
            break;
        }
    return ( estimate_plug() + (n-1.0)/n*f1 );
}


double Support::estimate_CL1() const // for small coefficient of variation
{
    double N1 = estimate_TG();
    double C = coverage_TG();

    double sum = 0.0;
    for ( const auto & pair : fin )
        sum += pair.first*( pair.first-1 )*pair.second;
    double gamma_sq = std::max( estimate_TG()*sum/n/(n-1) - 1 , 0.0 );
    return N1 + n*(1-C)/C*gamma_sq;
}


double Support::estimate_CL2() const // for large coefficient of variation
{
    double N1 = estimate_TG();
    double C = coverage_TG();
    
    double sum = 0.0;
    for ( const auto & pair : fin )
        sum += pair.first*( pair.first-1 )*pair.second;
    double gamma_sq = std::max( estimate_TG()*sum/n/(n-1) - 1 , 0.0 );
    double gamma_sq2 = std::max( gamma_sq*( 1 + n*(1-C)*sum/n/(n-1)/C ), 0.0 );
	
    return N1 + n*(1-C)/C*gamma_sq2;
}

/* ----------------------END OTHER ESTIMATORS-------------------------------- */



/* ----------------------OLD ESTIMATORS-------------------------------- */

double Support::estimate_old() const
{
    if ( pmin < Ratio / n ) // const
    {
        double result = 0.0;
        for ( const auto & pair : fin )
            if ( pair.first < N_thr ) 
            {
                result += getCoeff_old(pair.first) * pair.second;
            }
            else
            {
                result += pair.second;
            }
        return result;
    }
    else
        return estimate_plug();
}

double Support::getCoeff_old( double N ) const
{
    double A = ( Ratio + n*pmin ) / ( Ratio - n*pmin );
    double a[L+1];
	
    ChebMore cheb(L, 1, -A); // polynomial of cos L arccos(t-A)
    std::vector<double> a0 = cheb.expand(); // Expand: cos L arccos(x-A)=sum_i a0[i]*x^i
    double amp = cheb.evaluate(0); // cos L arccos(-A)

    // Expand: 1- cos L arccos(t-A)/cos L arccos(-A) = sum_i a[i]*x^i
    for (int i = 0;i < L+1;i++) a[i] = - a0[i] / amp;
    a[0] = 0;
    
    double s = 2 / ( Ratio - n * pmin );
    double gL = a[L];
    for (int i = L - 1; i>=0; i--)
        gL = a[i] + gL * (N-i) * s;
    return gL;
}



double Support::estimate_sinc() const
{
    double sum = 0.0;
    double k = 1.0/pmin;
    for ( const auto & pair : fin )
    {
        int j = pair.first;
        int fj = pair.second;
        if ( pair.first <= L )
        {
            double gj = 0;
            double multi = -1.0;
            for ( int i = 1; 2*i <= j; i++ )
            {
                multi *= - (3.1415926*k/n) * (3.1415926*k/n);
                gj += binom( j , 2*i ) * multi / (2*i+1);
            }
            sum += gj * fj;
        }
        else
        {
            sum += fj;
        }
    }
	
    return sum;
}

/* ----------------------END OLD ESTIMATORS-------------------------------- */



Support::Support()
{
}


Support::Support(double p_min)
    :Support()
{
    pmin = p_min;
    setDegree( 0.45*log(1.0/pmin) );      // Default L = 0.56 * log(k), where k = 1/pmin
    setInterval( 0.5*log(1.0/pmin) );     // Default interval [1/k, 0.5*log(k)/n], where k = 1/pmin
}


void Support::setFin(std::shared_ptr< const std::map<int, int> > ptr_fin_map) 
{
    fin.clear();
    n = 0;
    for ( auto it = ptr_fin_map->begin(); it != ptr_fin_map->end(); ++it )
    {
        int freq = it->first, cnt = it->second;
        fin.push_back( std::make_pair( freq, cnt ) );
        n += (freq * cnt);
    }
}

void Support::setFin(const std::string filename) 
{
    fin.clear();
    n = 0;
    std::ifstream infile;
    infile.open( filename.c_str() );
    int freq, cnt;
    while ( (infile>>freq).good() )
    {
        infile >> cnt;
        fin.push_back( std::make_pair( freq, cnt ) );
        n += (freq * cnt);
    }
    infile.close();
}

void Support::setFin(const std::vector<int> &freq_in, const std::vector<int> &cnt_in)
{
    fin.clear();
    n = 0;
    for ( int i = 0; i < freq_in.size(); i++)
    {
        int freq = freq_in[i], cnt = cnt_in[i];
        fin.push_back( std::make_pair( freq, cnt ) );
        n += (freq * cnt);
    }
}


void Support::setHist( const std::vector<unsigned> &hist )
{
    fin.clear();
    n = 0;
    
    std::map<int, int> fin_map;
    for ( const auto & freq : hist )
    {
        auto iter = fin_map.find( freq );
        if ( iter == fin_map.end() )
            fin_map.insert( std::make_pair( freq,1 ) );
        else
            ++(iter->second);
    }
    
    for ( auto it = fin_map.begin(); it != fin_map.end(); ++it )
    {
        int freq = it->first, cnt = it->second;
        fin.push_back( std::make_pair( freq, cnt ) );
        n += (freq * cnt);
    }
}

void Support::setHist(const std::string filename) 
{
    fin.clear();
    n = 0;
    std::map<int, int> fin_map;
    
    std::ifstream infile;
    infile.open( filename.c_str() );
    int freq;
    while ( (infile>>freq).good() )
    {
        if ( freq == 0)
            continue;

        auto iter = fin_map.find( freq );
        if ( iter == fin_map.end() )
            fin_map.insert( std::make_pair( freq,1 ) );
        else
            ++(iter->second);
    }
    infile.close();

    for ( auto it = fin_map.begin(); it != fin_map.end(); ++it )
    {
        int freq = it->first, cnt = it->second;
        fin.push_back( std::make_pair( freq, cnt ) );
        n += (freq * cnt);
    }
}
