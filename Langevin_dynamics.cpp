/*
Library of functions for lattice gas simulations and otherwise. C++ implemntations, making them classes. 
*/
//#define CLUSTER
#define PI 3.14159
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <sstream>
#include "Langevin_dynamics.h"
#ifdef CLUSTER
#include "/home/svaikunt/local/lib/include/fftw3.h"
#endif 
#ifndef CLUSTER
#include <fftw3.h>
#endif
#include <assert.h>
using namespace std;

Langevin_dynamics::~Langevin_dynamics(){
	gsl_rng_free(r);
}

Langevin_dynamics::Langevin_dynamics(){}



void Langevin_dynamics::initialize(int N_max1,long int randomseed1,double S1){
    Lx=10; //size of box
    Ly=10; //size of box
    N_type=2;// number of particle types
    Ntype[0]=4; //There are two particle of type 0
    Ntype[1]=N_max1;//Number of particles of type 1 input by user.
	S=S1; //value of S for biasing.
	N_max=N_max1;
    randomseed=randomseed1;
	r=gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(r,randomseed);
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
        fpos.push_back(fpos2);
        fpostype.push_back(fpos2);
        upostype.push_back(fpos2);
        zerovec.push_back(fpos2);
    }
    equilibrate();//initial equilibration run.
    //pos[type][Number][dimension];
}

void Langevin_dynamics::equilibrate(){
    double fd_term,noise_0,noise_1;
    double del1,del2;
    gamma_i=1;
    //Compute forces, evolve dynamics, return y factor
    double time=0;
    double dt=0.0001;
    int count;
    for (time=0;time<dt*800000;time=time+dt){
        compute_forces();
        count+=1;
        for (int loopi=0;loopi<N_type;loopi++){
            for (int loopj=0;loopj<Ntype[loopi];loopj++){
                fd_term = sqrt( 2 * dt / (gamma_i));
                noise_0 = fd_term*gsl_ran_gaussian(r,1);
                noise_1 = fd_term * gsl_ran_gaussian(r,1);
                del1=dt * fpos[loopi][loopj][0] + noise_0;
                del2=dt * fpos[loopi][loopj][1] + noise_1;
                if (fabs(del1)>0.25| fabs(del2)>0.25|| pos[loopi][loopj][0]>Lx || pos[loopi][loopj][0]<0||pos[loopi][loopj][1]>Ly||pos[loopi][loopj][1]<0){
                    cout<<del1<<"\t"<<noise_0<<"\t"<<dt * fpos[loopi][loopj][0]<<"Kill program\n";
                    cout<<del2<<"\t"<<noise_1<<"\t"<<dt * fpos[loopi][loopj][1]<<"Kill program\n";
                    cout<<"Time"<<"\t"<<time<<"\n";
                    //Check for big overlaps during equilibration and handle them.
                    if (del1>0.25)
                        del1=0.25;
                    if (del1<-0.25)
                        del1=-0.25;
                    if (del2>0.25)
                        del2=0.25;
                    if (del2<-0.25)
                        del2=-0.25;

                    pos[loopi][loopj][0] += del1;
                    pos[loopi][loopj][1] += del2;
                }
                else{
                    pos[loopi][loopj][0] += del1;
                    pos[loopi][loopj][1] += del2;
                }
                if (pos[loopi][loopj][0]>Lx)
                    pos[loopi][loopj][0]=pos[loopi][loopj][0]-Lx;
                if (pos[loopi][loopj][0]<0)
                    pos[loopi][loopj][0]+=Lx;
                if (pos[loopi][loopj][1]>Ly)
                    pos[loopi][loopj][1]=pos[loopi][loopj][1]-Ly;
                if (pos[loopi][loopj][1]<0)
                    pos[loopi][loopj][1]+=Ly;
            }
        }
    }

}

int Langevin_dynamics::propogate_dynamics(double dt){
    double fd_term,noise_0,noise_1;
    double del1,del2;
    gamma_i=1;
    //Compute forces, evolve dynamics, return y factor
    compute_forces();
    for (int loopi=0;loopi<N_type;loopi++){
        for (int loopj=0;loopj<Ntype[loopi];loopj++){
            fd_term = sqrt( 2 * dt / (gamma_i));
            noise_0 = fd_term*gsl_ran_gaussian(r,1);
            noise_1 = fd_term * gsl_ran_gaussian(r,1);
            del1=dt * fpos[loopi][loopj][0] + noise_0;
            del2=dt * fpos[loopi][loopj][1] + noise_1;
            if (fabs(del1)>0.5*Lx|| fabs(del2)>0.5*Lx){
                cout<<del1<<"\t"<<del2<<"Kill program\n";
                assert(0);//To check for wildly inappropriate moves.
            }
            pos[loopi][loopj][0] += del1;
            pos[loopi][loopj][1] += del2;
            if (pos[loopi][loopj][0]>Lx)
                pos[loopi][loopj][0]=pos[loopi][loopj][0]-Lx;
            if (pos[loopi][loopj][0]<0)
                pos[loopi][loopj][0]+=Lx;
            if (pos[loopi][loopj][1]>Ly)
                pos[loopi][loopj][1]=pos[loopi][loopj][1]-Ly;
            if (pos[loopi][loopj][1]<0)
                pos[loopi][loopj][1]+=Ly;
        }
    }
    double y1=computey();
    double y2=exp(-dt*S*y1);
    //cout<<y<<"\n";
    int yc;
    if (fabs(y1*S*dt)>5){
        cout<<"Panic in system"<<"\t"<<y1<<"\n";
        cout.flush();
        return 0.0;
        //Absurd value of y are not returned.
    }
    else{
        yc=(int)(floor(y2+gsl_rng_uniform(r)));
        cout.flush();
        return yc;//yc is used to compute cloning probabilities. 
    }
}

void Langevin_dynamics::compute_forces(){
    double f12x=0;
    double f12y=0;
    double x12=0;
    double y12=0;
    double r12=0;
    double f122x=0;
    double f122y=0;
    fpos=zerovec;
    fpostype=zerovec;
    upostype=zerovec;
    for (int loopi=0;loopi<N_type;loopi++){
        for(int loopk=0;loopk<Ntype[loopi];loopk++){
            for(int loopj=0;loopj<N_type;loopj++){
                for(int loopl=0;loopl<Ntype[loopj];loopl++){
                    f12x=0;
                    f12y=0;
                    f122x=0;
                    f122y=0;
                    if (loopi!=loopj || loopk!=loopl){
                        x12=pos[loopi][loopk][0]-pos[loopj][loopl][0];
                        y12=pos[loopi][loopk][1]-pos[loopj][loopl][1];
                        if (x12>0.5*Lx)
                            x12=Lx-x12;
                        if (x12<-0.5*Lx)
                            x12=x12+Lx;
                        if (y12>0.5*Ly)
                            y12=Ly-y12;
                        if (y12<-0.5*Ly)
                            y12=y12+Ly;
                        r12=pow(x12*x12+y12*y12,0.5);
                        f12x=4*(12*pow(r12,-14)-6*pow(r12,-8))*x12;
                        f12y=4*(12*pow(r12,-14)-6*pow(r12,-8))*y12;
                        f122x=-4*(12*pow(r12,-13)-6*pow(r12,-7))*pow(r12,-1.0)+4*(168*pow(r12,-16)-48*pow(r12,-10.0))*(x12)*x12;
                        f122y=-4*(12*pow(r12,-13)-6*pow(r12,-7))*pow(r12,-1.0)+4*(168*pow(r12,-16)-48*pow(r12,-10.0))*(y12)*y12;
                        if (r12!=0 && r12<pow(2,1.0/6.0)){
                         fpos[loopi][loopk][0]+=f12x;
                         fpos[loopi][loopk][1]+=f12y;
                            if(r12<0.75){
                                cout<<"Panic\t"<<pos[loopi][loopk][0]<<"\t"<<pos[loopj][loopl][0]<<"\t"<<pos[loopi][loopk][1]<<"\t"<<pos[loopj][loopl][1]<<"\t"<<f12x<<"\t"<<f12y<<"\n";
                                //Just to check for big overlaps. Might happen during equilibration.
                            }
                        }
                        if (loopi!=loopj){
                            if (r12!=0 && r12<pow(2,1.0/6.0)){
                             fpostype[loopi][loopk][0]+=f12x;
                             fpostype[loopi][loopk][1]+=f12y;
                             upostype[loopi][loopk][0]+=f122x;
                             upostype[loopi][loopk][1]+=f122y;
                            }
                        }
                    }
                }
            }
        }
    }
}

double Langevin_dynamics::computey(){
    double y[2];
    y[0]=0;
    y[1]=0;
    double xy[2];
    xy[0]=0;
    xy[1]=0;
    for (int loopi=0;loopi<N_type;loopi++){
        for(int loopj=0;loopj<Ntype[loopi];loopj++){
            xy[loopi]+=fpos[loopi][loopj][0]*(-fpostype[loopi][loopj][0])+fpos[loopi][loopj][1]*(-fpostype[loopi][loopj][1]);
            y[loopi]+=fpos[loopi][loopj][0]*(-fpostype[loopi][loopj][0])+fpos[loopi][loopj][1]*(-fpostype[loopi][loopj][1])+upostype[loopi][loopj][0]+upostype[loopi][loopj][1];
        }
    }
    return y[0]+y[1];//yfactor computed for cloning.
}



