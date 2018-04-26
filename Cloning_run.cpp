
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
    double y[N_snapshots];
    double yc[N_snapshots];
    double yc2[N_snapshots];
    double sumy=0;
	vector< vector < vector < vector<double> > > > tpos;
    
    for (int i=0;i<N_snapshots;i++){
        vector< vector< vector<double> > > pos;
        for (int loopi=0;loopi<N_type;loopi++){
            vector< vector<double> > pos2;
            vector< vector<double> >fpos2;
            for (int loopj=0;loopj<Ntype[loopi];loopj++){
                vector<double> pos1;
                vector<double> fpos1;
                for(int loopk=0;loopk<2;loopk++){
                    pos1.push_back(gsl_rng_uniform(r)*Lx);
                    fpos1.push_back(0);
                }
                pos2.push_back(pos1);
                fpos2.push_back(fpos1);
            }
            pos.push_back(pos2);
        }
        tpos.push_back(pos);
    }
    
    for (int i=0;i<N_snapshots;i++){
        cout<<"Clone>>>>>>>>>>>>>>>\t"<<i<<"\n";
		mylattice[i].initialize(N_max,randomseed+i,S);
        tpos[i]=mylattice[i].pos;
    }
    
	for (double teetotaler=0;teetotaler<t_analysis;teetotaler+=dt){
        sumy=0;
        double sumyc=0;
        for (int loopj=0;loopj<N_snapshots;loopj++){
            //j=pick randomly from the clones
            y[loopj]=mylattice[loopj].propogate_dynamics(dt);//propogate clone.
            tpos[i]=mylattice[i].pos;//storing clone positions for switching.
            sumy+=y[loopj];
        }
        for (int loopj=0;loopj<N_snapshots;loopj++){
            yc[loopj]=(int)(floor(y[loopj]/sumy*N_snapshots+gsl_rng_uniform(r)));
            sumyc+=yc[loopj];
            yc2[loopj]=sumyc;
        }
        for (int loopj=0;loopj<N_snapshots;loopj++){
            double randtemp=gsl_rng_uniform(r)*sumyc;
            int iterate=0;
            do{
                iterate+=1;
            }while(yc[iterate]<=randtemp);
            mylattice[loopj].pos=tpos[iterate];
        }
        double ratio=((double)(sumy)*pow(N_snapshots,-1.0);
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
        fileout<<N_max+4<<"\n";
        fileout<<"Next\n";
        for (int i=0;i<N_max;i++){
        fileout<<"H"<<"\t"<<mylattice[loopi].pos[1][i][0]<<"\t"<<mylattice[loopi].pos[1][i][1]<<"\t"<<0.1<<"\n";
        }
        fileout<<"O"<<"\t"<<mylattice[loopi].pos[0][0][0]<<"\t"<<mylattice[loopi].pos[0][0][1]<<"\t"<<0.1<<"\n";
        fileout<<"O"<<"\t"<<mylattice[loopi].pos[0][1][0]<<"\t"<<mylattice[loopi].pos[0][1][1]<<"\t"<<0.1<<"\n";
        fileout<<"O"<<"\t"<<mylattice[loopi].pos[0][2][0]<<"\t"<<mylattice[loopi].pos[0][2][1]<<"\t"<<0.1<<"\n";
        fileout<<"O"<<"\t"<<mylattice[loopi].pos[0][3][0]<<"\t"<<mylattice[loopi].pos[0][3][1]<<"\t"<<0.1<<"\n";
        
    }
   
	
	
}

