#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "support.h"
#include "commandline.h"

void print_param(const Support &support);
void print_results(const Support &support);


int main(int argc, char *argv[])
{
    double pmin = 0;
    int L = 0;
    double M = 0;
    std::string fin = "", hist = "";
    
    std::CommandLine cmd;
    cmd.AddValue ("pmin",  "Minimum mass", pmin);
    cmd.AddValue ("L",  "Polynoimal degree. Default c0*log(1/pmin)", L);
    cmd.AddValue ("M",  "M/n is the right-end of approximation interval. Default c1*log(1/pmin)", M);
    cmd.AddValue ("fin",  "fingerprint data file", fin);
    cmd.AddValue ("hist",  "histogram data file", hist);
    cmd.Parse (argc, argv);
    if( pmin == 0 ) { std::cerr<<"Please input pmin!\n"; exit(1); }
    if( L == 0 ) L = 0.45*log(1/pmin); // Default value
    if( M == 0 ) M = 0.5*log(1/pmin);  // Default value
    if( (hist=="") && (fin=="") ) { std::cerr<<"Please input fingerprint or histogram!\n"; exit(1); }
    
    printf("\n");
    // Configure estimator
    Support support( pmin );  // set pmin
    support.setDegree( L );   // set degree    =L
    support.setInterval( M ); // set interval  =[pmin,M/n]
    print_param(support);

    // Input fingerprint from file
    if ( fin!="")
        support.setFin(fin);
    else
        support.setHist(hist);
    print_results(support);

    
    return 0;
}


void print_param(const Support &support)
{
    printf("Parameters:\n");
    printf("Preset pmin value\t=%.2e\n", support.getPmin());
    printf("Degree of polynomial\t=%d\n", support.getDegree());
    printf("Approximation interval\t=[%.2e,%.2f/n]\n", support.getPmin(), support.getInterval());
    printf("\n");
}


void print_results(const Support &support)
{
    printf("Results:\n");
    printf("Sample size\t=%d\n",support.getSampleSize());
    printf("Polynomial\t=%d \n",(int)support.estimate());
    printf("Plug-in\t\t=%d \n",(int)support.estimate_plug());
    printf("\n");
}


