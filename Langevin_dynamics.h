/*
Library of functions for lattice KMC simulations and otherwise. C++ implemntations, making them classes. 
*/
#ifndef SOSGAURD
#define SOSGAURD
#include <cstdlib>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <fstream>
#include <vector>

using namespace std;


class Langevin_dynamics{
private:
	double t_analysis;
	double Lx;
	double Ly;
	double Ntype[2];
	double N_type;
	int N_max;
	double S;
	double gamma_i;
public:
  
  ofstream fileout_traj,fileout_activity;
  
  
  gsl_rng *r;
  long int randomseed;
  
  vector < vector < vector<double> > > pos;
  vector < vector < vector<double> > > fpos;
  vector < vector < vector<double> > > fpostype; 
  vector < vector < vector<double> > > upostype;
  vector < vector < vector<double> > > zerovec;

  
  Langevin_dynamics();
  ~Langevin_dynamics();
  void initialize(int, long int,double);
  double propogate_dynamics(double);
  void compute_forces();
  double computey();
  void equilibrate();
};

#endif



