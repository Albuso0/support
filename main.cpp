#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "support.h"
#include "mathmore.h"
#include "commandline.h"


int main(int argc, char *argv[])
{
    double pmin = 1.0/1e5;
    double cL = 0.56, cp = 0.5;
    std::CommandLine cmd;
    cmd.AddValue ("pmin",  "pmin", pmin);
    cmd.AddValue ("cL",  "L=cL log k", cL);
    cmd.AddValue ("cp",  "rEnd=cp log k/n", cp);
    cmd.Parse (argc, argv);

    Support support( pmin ); // set pmin
    support.setInterval( cL*log(1.0/pmin) ); // Approximation interval is [1/k, interval/n ]
    support.setDegree( cp*log(1.0/pmin)); // Polynomial degree. Plug-in if N>L


    support.setFin( "fin_sample.txt" );
    std::cout<<support.estimate()<<std::endl;

    
    return 0;
}


