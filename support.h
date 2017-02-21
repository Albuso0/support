#ifndef SUPPORT_H
#define SUPPORT_H

#include <map>
#include <vector>

class Support
{
public:
    Support();
    Support(double pmin);
    virtual ~Support(){}

	
    // Direct design, no sample splitting in theory
    double estimate() const;
    double getCoeff( int N ) const;
	
    // Other estimators
    double estimate_plug() const;  // Plug-in estimator
    double coverage_TG() const;    // Turing-Good coverage
    double estimate_TG() const;    // Turing-Good estimator
    double estimate_J1() const;    // First order Jackknife
    double estimate_Chao1() const;   // Chao 1
    double estimate_CL1() const;   // Chao-Lee 1
    double estimate_CL2() const;   // Chao-Lee 2

	
    void setPmin( double p_min ) { pmin = p_min; }
    void setInterval( double rEnd ){ Ratio = rEnd; }
    void setDegree( int deg ) { L = deg; }
    
    // set fingerprint, also update sample size
    void setFin( const std::string filename );
    void setFin( const std::vector<int> &freq, const std::vector<int> &cnt );
    // set fingerprint through histogram, also update sample size
    void setHist( const std::string filename );
    void setHist( const std::vector<unsigned> &hist );
    

    double getPmin() const{ return pmin; }
    int getDegree() const{ return L; }
    double getInterval() const{ return Ratio; }
    int getSampleSize() const{ return n; }
private:
    double pmin;  // =1/k. Preset minimum non-zero mass
    int L;        // =c0*log(k). Degree of polynomial
    double Ratio; // =c1*log(k). Approximation over [pmin,c1*log(k)/n]
    int n;        // sample size
    	
    std::vector< std::pair<int, int> > fin; // Fingerprint(profile). Sorted fingerprint preferred
    
};


#endif

