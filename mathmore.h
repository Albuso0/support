#ifndef MATHMORE_H
#define MATHMORE_H
#include <memory>
#include <vector>


long long binom(unsigned n, unsigned k);

double rmse( const std::vector<int> &v, int truth );
double stdev( const std::vector<int> &v );
double mean( const std::vector<int> &v );
double min_positive_normalized( const std::vector<double> &p ); // 0<=output<=1
double cnt_positive( const std::vector<double> &p );

class ChebMore  // Polynomial in the form cos(L arccos (bx+c))=sum a_i * x^i
{
public:
    ChebMore(){}
    ChebMore(int _L, double _b = 1, double _c = 0):L(_L),b(_b),c(_c){}
    virtual ~ChebMore(){}
    

    std::vector<double> expand_basic() const; // Expansion of cos(L arccos x)
    std::vector<double> expand() const;       // Expansion of cos(L arccos (bx+c))
    
    double evaluate(double x) const; // Evaluation of cos(L arccos (bx+c))

private:
    int L;
    double b,c;
};

#endif
