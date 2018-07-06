
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
	  cout<<"Usage N_max(number of particles of type1, there are 4 of type 0) t_analysis randomseed timestep  S(biasing function) Pe tau soft ksoft_solute_solvent N_Snapshots\n";
	  exit(1);
	}
	int N_max=atoi(argv[1]);
	double t_analysis=atof(argv[2]);
	long int randomseed=atoi(argv[3]);
	double dt=atof(argv[4]);
	double S=atof(argv[5]);
	double Pe=atof(argv[6]);
	double tau=atof(argv[7]);
	int soft=atoi(argv[8]);
	double ksoft=atof(argv[9]);
	int N_snapshots=atoi(argv[10]);//number of clones
        double Lx=7.0;
        double Ly=7.0;
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
    
    
    int N_type=2;// number of particle types
    int Ntype[2];
    Ntype[0]=4; //There are 4  particle of type 0
    Ntype[1]=N_max;//Number of particles of type 1 input by user.
    for (int i=0;i<N_snapshots;i++){
        vector< vector< vector<double> > > pos;
        for (int loopi=0;loopi<N_type;loopi++){
            vector< vector<double> > pos2;
            for (int loopj=0;loopj<Ntype[loopi];loopj++){
                vector<double> pos1;
                for(int loopk=0;loopk<2;loopk++){
                    pos1.push_back(0);
                }
                pos2.push_back(pos1);
            }
            pos.push_back(pos2);
        }
        tpos.push_back(pos);
    }
    
    Langevin_dynamics masterlattice;
    masterlattice.initialize(N_max,Ntype[0],randomseed,S,soft,ksoft);
    cout<<"Clone>>>>>>>>>>>>>>>Equilibrating master\n";
    masterlattice.equilibrate(Pe,tau);
    
    for (int i=0;i<N_snapshots;i++){
        cout<<"Clone>>>>>>>>>>>>>>>\t"<<i<<"\n";
        for (int j=0;j<5.0/dt;j++)
          double dump=masterlattice.propogate_dynamics(dt,Pe,0,tau);
        mylattice[i].initialize(N_max,Ntype[0],randomseed+i,S,soft,ksoft);
        //mylattice[i].equilibrate(Pe,tau);
         mylattice[i].pos=masterlattice.pos;
        tpos[i]=mylattice[i].pos;
    }
    
    char outputfile1[100];
    sprintf(outputfile1,"g_rN%d.N2%d.S%.3f.Lx%.2f.Clone%d.Soft%d.ktracer%.2f.E10.XYZ",N_max,Ntype[0],S,Lx,N_snapshots,soft,ksoft);
    ofstream fileoutg;
    fileoutg.open(outputfile1);
    double xdis,ydis,r2;
    int countdt=0;
    int dumpdt=(int)(5.0/dt);
	for (double teetotaler=0;teetotaler<t_analysis;teetotaler+=dt){
        countdt+=1;
        sumy=0;
        double sumyc=0;
        if (countdt%dumpdt==0 && countdt>dumpdt){
            for (int loopi=0;loopi<N_snapshots;loopi++){
                for( int i=0; i<Ntype[0];i++){
                    for(int j=0;j<N_max;j++){
                        xdis=mylattice[loopi].pos[0][i][0]-mylattice[loopi].pos[1][j][0];
                        ydis=mylattice[loopi].pos[0][i][1]-mylattice[loopi].pos[1][j][1];
                        if (xdis>0.5*Lx)
                            xdis=xdis-Lx;
                        if (xdis<-0.5*Lx)
                            xdis=xdis+Lx;
                        if (ydis>0.5*Ly)
                            ydis=ydis-Ly;
                        if (ydis<-0.5*Ly)
                            ydis=ydis+Ly;
                        r2=pow(xdis*xdis+ydis*ydis,0.5);
                        fileoutg<<r2<<"\n";
                    }
                }
            }
        }
        for (int loopj=0;loopj<N_snapshots;loopj++){
            //j=pick randomly from the clones
            y[loopj]=mylattice[loopj].propogate_dynamics(dt,Pe,teetotaler,tau);//propogate clone.
            
            tpos[loopj]=mylattice[loopj].pos;//storing clone positions for switching.
            sumy+=y[loopj];
        }
        for (int loopj=0;loopj<N_snapshots;loopj++){
            yc[loopj]=(int)(floor(y[loopj]/sumy*N_snapshots+gsl_rng_uniform(r)));
            //cout<<"Testing yc\t"<<yc[loopj]<<"\n";
            sumyc+=yc[loopj];
            yc2[loopj]=sumyc;
        }
        if(sumyc==N_snapshots){
            int clonecount=0;
            for (int loopj=0;loopj<N_snapshots;loopj++){
                for (int loopk=0;loopk<yc[loopj];loopk++){
                    mylattice[loopk+clonecount].pos=tpos[loopj];
                }
                clonecount+=yc[loopj];
                //cout<<loopj<<"\t"<<clonecount<<"\n";
                //cout.flush();
            }
        }
        if (sumyc!=N_snapshots){
            for (int loopj=0;loopj<N_snapshots;loopj++){
                double randtemp=gsl_rng_uniform(r)*sumyc;
                int iterate=-1;
                do{
                    iterate+=1;
                }while(yc2[iterate]<=randtemp);
                if (iterate>N_snapshots-1){
                    cout<<"Iterate\t"<<iterate<<"\t"<<randtemp<<"\t"<<yc[N_snapshots-1]<<"\n";
                    iterate=N_snapshots-1;
                }
                mylattice[loopj].pos=tpos[iterate];
                if(teetotaler>0.0*t_analysis)
                    cout<<"Clone\t"<<loopj<<"\t is now clone \t"<<iterate<<"\t"<<yc2[iterate]<<"\t"<<randtemp<<"\n";
            }
        }
        double ratio=((double)(sumy)*pow(N_snapshots,-1.0));
        growthcgf=growthcgf*ratio;
        //compute y for clone. Kill (copy some other clone into this) or select clones randomly to copy into
        //keep track of growth function.
    }
    cout<<"CGF:"<<log(growthcgf)/t_analysis<<"\n";
    
    //Compute averages over the cloned lattices.
    double avgy=0; //Average value of dU12/dt will be stored in avgy.
    char outputfile[100];
    sprintf(outputfile,"ystatsNmax%d.N2%d.S%.3f.Lx%.2f.Clone%d.Soft%d.E10.ktracer%.2f.XYZ",N_max,Ntype[0],S,Lx,N_snapshots,soft,ksoft);
    ofstream fileoutystats;
    fileoutystats.open(outputfile);
    fileoutystats<<"CGF\t"<<log(growthcgf)/t_analysis<<"\n";
    for (int loopi=0;loopi<N_snapshots;loopi++){
        mylattice[loopi].compute_forces();
        double tempy=mylattice[loopi].computey(Pe,t_analysis,tau);
        fileoutystats<<tempy<<"\n";
        cout<<"Compute y is "<<tempy<<"\t Clone is \t"<<loopi<<"\n";
        avgy+=tempy;
    }
	cout<<"avgenergy:"<<avgy*pow(N_snapshots,-1.0)<<"\n";
    sprintf(outputfile,"SnapshotsN%d.N2%d.S%.3f.Lx%.2f.Clone%d.Soft%d.E10.ktracer%.2f.XYZ",N_max,Ntype[0],S,Lx,N_snapshots,soft,ksoft);
    ofstream fileout;
    fileout.open(outputfile);
    for (int loopi=0;loopi<N_snapshots;loopi++){
        fileout<<N_max+Ntype[0]<<"\n";
        fileout<<"Next\n";
        for (int i=0;i<N_max;i++){
            fileout<<"H"<<"\t"<<mylattice[loopi].pos[1][i][0]<<"\t"<<mylattice[loopi].pos[1][i][1]<<"\t"<<0.1<<"\n";
        }
        for (int i=0;i<Ntype[0];i++){
            fileout<<"O"<<"\t"<<mylattice[loopi].pos[0][i][0]<<"\t"<<mylattice[loopi].pos[0][i][1]<<"\t"<<0.1<<"\n";
        }
    }
    
    
}

