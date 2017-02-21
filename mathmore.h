#ifndef MATHMORE_H
#define MATHMORE_H

#include <vector>

long long binom(unsigned n, unsigned k);

class ChebMore  // Polynomial in the form cos(L arccos (bx+c))=sum a_i * x^i
{
public:
    ChebMore(){}
    ChebMore(int L_in, double b_in = 1, double c_in = 0):L(L_in),b(b_in),c(c_in){}
    virtual ~ChebMore(){}
    

    std::vector<double> expand_basic() const; // Expansion of cos(L arccos x)
    std::vector<double> expand() const;       // Expansion of cos(L arccos (bx+c))
    
    double evaluate(double x) const; // Evaluation of cos(L arccos (bx+c))

private:
    int L;
    double b,c;
};

#endif
