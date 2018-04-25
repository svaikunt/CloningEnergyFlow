
/*
KMC simulation. */

#define PI 3.14159
//#define CLUSTER
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#ifdef CLUSTER
#include "/home/svaikunt/local/lib/include/fftw3.h"
#endif 
#ifndef CLUSTER
#include <fftw3.h>
#endif
#include <assert.h>
#include "Lattice_KMC.h"

#define ARGS 8

using namespace std;
int main( int argc,char *argv[]){

  int loopi,loopj;
  
	if (argc==1){
	  cout<<"KMC simulation\n";
	  cout<<"Usage N_max t_analysis randomseed N_snapshots S\n";
	  exit(1);
	}
	int N_max=atoi(argv[1]);
	double t_analysis=atof(argv[2]);
	long int randomseed=atoi(argv[3]);
	int N_snapshots=atoi(argv[4]);
	double S=atof(argv[5]);
	Lattice_KMC mylattice;
	mylattice.initialize(N_max,t_analysis,randomseed,S);
	mylattice.propagate(N_snapshots);
	
	
	
}

