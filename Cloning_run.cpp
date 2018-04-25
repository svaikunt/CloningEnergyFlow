
/*
CLONING simulation. */

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
#include "Langevin_dynamics.h"

#define ARGS 8

using namespace std;
int main( int argc,char *argv[]){

    int loopi,loopj;
    double growthcgf=1;;
    if (argc==1){
	  cout<<"Cloning simulation\n";
	  cout<<"Usage N_max(number of particles) t_analysis randomseed timestep  S(biasing function)\n";
	  exit(1);
	}
	int N_max=atoi(argv[1]);
	double t_analysis=atof(argv[2]);
	long int randomseed=atoi(argv[3]);
	double dt=atof(argv[4]);
	double S=atof(argv[5]);
	int N_snapshots=50;
    
    //Initiating gsl
    gsl_rng *r=gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(r,randomseed);
    

    
	//Iniatiate N_snapshots or clones of mylattice.
	Langevin_dynamics mylattice[N_snapshots];
	vector < vector < vector<double> > > tpos;
    
    for (int i=0;i<N_snapshots;i++){
        cout<<"Clone>>>>>>>>>>>>>>>\t"<<i<<"\n";
		mylattice[i].initialize(N_max,randomseed+i,S);
    }
    
	for (double teetotaler=0;teetotaler<t_analysis;teetotaler+=dt){
        for (int loopj=0;loopj<N_snapshots;loopj++){
            int j=gsl_rng_uniform_int(r,N_snapshots);
            //j=pick randomly from the clones
            int y=mylattice[j].propogate_dynamics(dt);//propogate clone.
                if (y==0){
                    do{
                        loopi=gsl_rng_uniform_int(r,N_snapshots);
                    }
                    while(loopi==j);
                    mylattice[j].pos=mylattice[loopi].pos;
                }
                else{
                    tpos=mylattice[j].pos;
                        for(loopi=0;loopi<y-1;loopi++){
                            loopj=gsl_rng_uniform_int(r,N_snapshots+y-1);
                            if (loopj<N_snapshots)
                                mylattice[loopj].pos=tpos;
                        }
                }
            double ratio=((double)(N_snapshots+y-1))*pow(N_snapshots,-1.0);
            growthcgf=growthcgf*ratio;
            //compute y for clone. Kill (copy some other clone into this) or select clones randomly to copy into
            //keep track of growth function.
        }
    }
    cout<<"CGF:"<<log(growthcgf)/t_analysis<<"\n";
    
    //Compute averages over the cloned lattices.
    double avgy=0; //Average value of dU12/dt will be stored in avgy.
    char outputfile[100];
    sprintf(outputfile,"ystatsNmax%d.S%.2f.XYZ",N_max,S);
    ofstream fileoutystats;
    fileoutystats.open(outputfile);
    fileoutystats<<"CGF\t"<<log(growthcgf)/t_analysis<<"\n";
    for (int loopi=0;loopi<N_snapshots;loopi++){
        mylattice[loopi].compute_forces();
        double tempy=mylattice[loopi].computey();
        fileoutystats<<tempy<<"\n";
        cout<<"Compute y is "<<tempy<<"\t Clone is \t"<<loopi<<"\n";
        avgy+=tempy;
    }
	cout<<"avgenergy:"<<avgy*pow(N_snapshots,-1.0)<<"\n";
    sprintf(outputfile,"SnapshotsN%d.S%.2f.XYZ",N_max,S);
    ofstream fileout;
    fileout.open(outputfile);
    for (int loopi=0;loopi<N_snapshots;loopi++){
        fileout<<N_max+2<<"\n";
        fileout<<"Next\n";
        for (int i=0;i<N_max;i++){
        fileout<<"H"<<"\t"<<mylattice[loopi].pos[1][i][0]<<"\t"<<mylattice[loopi].pos[1][i][1]<<"\t"<<0.1<<"\n";
        }
        fileout<<"O"<<"\t"<<mylattice[loopi].pos[0][0][0]<<"\t"<<mylattice[loopi].pos[0][0][1]<<"\t"<<0.1<<"\n";
        fileout<<"O"<<"\t"<<mylattice[loopi].pos[0][1][0]<<"\t"<<mylattice[loopi].pos[0][1][1]<<"\t"<<0.1<<"\n";
    }
   
	
	
}

