#include "support.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>

void print_param(const Support &support);
void print_results(const Support &support);

int main()
{
    double pmin = 1e-5;

    // Configure estimator
    Support support( pmin );
    // Optional parameters configuration
    support.setDegree( 5 ); 
    support.setInterval( 6 ); 
    print_param(support); // print parameters

    // Input fingerprint from file
    support.setFin("fin_sample.txt"); // fingerprint of 3000 samples generated from uniform[100000]
    print_results(support); // print estimation results

    // Input fingerprint inline
    std::vector<int> freq, cnt;
    freq.push_back(1); cnt.push_back(2910); 
    freq.push_back(2); cnt.push_back(45);
    support.setFin(freq, cnt);
    print_results(support);

    // Input histogram from file
    support.setHist("hist_sample.txt"); // histogram of 3000 samples generated from uniform[100000]
    print_results(support);
    
    // Input histogram inline
    std::vector<unsigned> hist;
    for (int i = 0; i < 2910; i++) hist.push_back(1);
    for (int i = 0; i < 45; i++) hist.push_back(2);
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


