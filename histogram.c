#include <time.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cj_utilities.c"

#define Err1 {printf("Please specify options correctly!\n"); exit(1);}

#define ARGS 4
#define LEN 1000

char *GetItemFromLine(char *line, int n);
int WordCount(char *line);
int CharacterType(char c);
int Bin(double x,double bsize);
void Usage(int argc, char **argv);
void JustTheStats(int argc, char **argv, int column);

main(argc,argv)
int argc;
char *argv[];
{
    int i,j,bmin,bmax,bin,nbin,column=1,PlotIt=0,ErrorBars=0,ct=1,argc0;
    double x,xmin,xmax,bsize,norm=0.,*count,*bars,avg=0.,std=0.;
    char line[LEN],*s,*Graphit,**argv0;
    FILE *fp_in,*fp_out;
	int skip=0;
	int skipcount=0;
    argc0=argc;
    argv0=malloc(argc*sizeof(char *));
    for (i=0;i<argc;i++) argv0[i] = argv[i];

    i=1;
    while (i<argc){
        if (strcmp(argv[i],"-c")==0){
            if (i==argc-1) Err1
            column = atoi(argv[i+1]);
            for (j=i+2;j<argc;j++) argv[j-2]=argv[j];
            argc -= 2;
        }
		else if (strcmp(argv[i],"-skip")==0){
            if (i==argc-1) Err1
				skip = atoi(argv[i+1]);
            for (j=i+2;j<argc;j++) argv[j-2]=argv[j];
            argc -= 2;
        }
		else if (strcmp(argv[i],"-n")==0){
            if (i==argc-1) Err1
            norm = atof(argv[i+1]);
            for (j=i+2;j<argc;j++) argv[j-2]=argv[j];
            argc -= 2;
        }
        else if (strcmp(argv[i],"-eb")==0){
            ErrorBars = 1;
            for (j=i+1;j<argc;j++) argv[j-1]=argv[j];
            argc -= 1;
        }
        else if (strcmp(argv[i],"!")==0){
            PlotIt = 1;
            for (j=i+1;j<argc;j++) argv[j-1]=argv[j];
            argc -= 1;
        }
        else i++;
    }

    if (argc==2) JustTheStats(argc,argv,column);

    Usage(argc,argv);

    fp_in  = Open(argv[1],"r");
    fp_out = Open(argv[2],"w");
    bsize  = atof(argv[3]);

    ReadComments(fp_in);

    fgets(line,LEN,fp_in);
    xmin = xmax = x = atof(GetItemFromLine(line,column-1));
    bmin = bmax = Bin(x,bsize);

    while (fgets(line,LEN,fp_in) != NULL){
        x = atof(GetItemFromLine(line,column-1));
		if (skipcount>=skip){
			ct += 1;
			bin = Bin(x,bsize);
			if (bin<bmin) bmin=bin;
			if (bin>bmax) bmax=bin;
			if (x<xmin) xmin=x;
			if (x>xmax) xmax=x;
		}
		skipcount=skipcount+1;
    }
	skipcount=0;
    fclose(fp_in);
    fp_in = Open(argv[1],"r");
    ReadComments(fp_in);

    nbin = bmax-bmin+3;
    count = malloc(nbin*sizeof(double));
    for (i=0;i<nbin;i++) count[i]=0.;

    bmin--;

    while (fgets(line,LEN,fp_in) != NULL){
        x = atof(GetItemFromLine(line,column-1));
		skipcount=skipcount+1;
		if(skipcount>=skip){
			count[Bin(x,bsize)-bmin]+=1.0;
			avg += x;
			std += x*x;
		}
    }

    avg /= ct;
    std = sqrt( (std/ct) - avg*avg );

    if (ErrorBars){
        bars = malloc(nbin*sizeof(double));
        for (i=0;i<nbin;i++) bars[i] = sqrt(count[i]);
    }

    printf("Input file:     %s\tcolumn %d\n",argv[1],column);
    printf("Output file:    %s\n",argv[2]);
    printf("Bin size:       %f\n",bsize);
    if (norm) printf("Normalization:  %f\n", norm);
    else      printf("Normalization:  %d\n", ct);
    printf("Average:        %f\n",avg);
    printf("Std. Deviation: %f\n",std);
    printf("Min,Max values: %f,%f\n",xmin,xmax);

    PrintCommandLine(fp_out,argc0,argv0);
    fprintf(fp_out,"#\n");
    fprintf(fp_out,"# Input file:     %s\tcolumn %d\n",argv[1],column);
    fprintf(fp_out,"# Output file:    %s\n",argv[2]);
    fprintf(fp_out,"# Bin size:       %f\n",bsize);
    if (norm) fprintf(fp_out,"# Normalization:  %f\n", norm);
    else      fprintf(fp_out,"# Normalization:  %d\n", ct);
    fprintf(fp_out,"# Average:        %f\n",avg);
    fprintf(fp_out,"# Std. Deviation: %f\n",std);
    fprintf(fp_out,"# Min,Max values: %f,%f\n#\n",xmin,xmax);
    fprintf(fp_out,"# NOTE: The first column gives the CENTER of the bin\n#\n");

    if (norm){
        for (i=0;i<nbin;i++) count[i] *= (norm/(bsize*ct));
        if (ErrorBars)
            for (i=0;i<nbin;i++) bars[i] *= (norm/(bsize*ct));
    }
    for (i=0;i<nbin;i++)
        if (ErrorBars)
            fprintf(fp_out,"%f\t%f\t%f\n",(bmin+i)*bsize+0.5*bsize,count[i],bars[i]);
        else
            fprintf(fp_out,"%f\t%f\n",(bmin+i)*bsize+0.5*bsize,count[i]);
    fclose(fp_out);

    if (PlotIt){
        Graphit = malloc((strlen(argv[2])+29)*sizeof(char));
        strcpy(Graphit,"python ~/Graphit/Graphit.py ");
        strcat(Graphit,argv[2]);
        system(Graphit);
    }

}

/*----------------------------------------------------------------*/

char *GetItemFromLine(char *line, int n)

{
    int i=0,m=0,len=0,p;
    char *ptr;

    if (WordCount(line)<n+1){
        printf("Error: small word count\n%s\n",line);
        exit(1);
    }

    while (CharacterType(line[i])==0) i++;
    p = i;

    while (m<n){
        while (CharacterType(line[i])==1) i++;
        while (CharacterType(line[i])==0) i++;
        p = i;
        m++;
    }

    while (CharacterType(line[p+len])==1) len++;

    ptr = malloc((len+1)*sizeof(char));
    for (i=0;i<len;i++) ptr[i]=line[p+i];
    ptr[len] = '\0';

    return ptr;
}

/*----------------------------------------------------------------*/

int WordCount(char *line)

{
    int wc=0,i=0;

    while (CharacterType(line[i])==0) i++;
    while (CharacterType(line[i])!=2){
        wc++;
        while (CharacterType(line[i])==1) i++;
        while (CharacterType(line[i])==0) i++;
    }

    return wc;
}

/*----------------------------------------------------------------*/

int CharacterType(char c)

{
    switch (c){
        case ' '  : return 0;
        case '\t' : return 0;
        case '\n' : return 2;
        case EOF  : return 2;
        default   : return 1;
    }
}

/*----------------------------------------------------------------*/

int Bin(double x,double bsize)

{
    if (x>0.) return (int) (x/bsize);
    else      return (int) (x/bsize)-1;
}

/*----------------------------------------------------------------*/

void Usage(int argc, char **argv)

{
        if(argc!=ARGS){
        printf("Usage: %s infile outfile binsize","Histogram");
        printf(" [-c column] [-skip skip] [-n normalization] [-eb (error bars)] [!]\n");
        exit(1);
    }
}

/*----------------------------------------------------------------*/

void JustTheStats(int argc, char **argv, int column)

{
    int ct,first=1;
    double x,xmin,xmax,avg=0.,std=0.;
    char line[LEN];
    FILE *fp_in;

    ReadComments( fp_in = Open(argv[1],"r") );

    while (fgets(line,LEN,fp_in) != NULL){
        x = atof(GetItemFromLine(line,column-1));
        if (first){ xmin=xmax=x; first=0; }
        else{ if (x<xmin) xmin=x; if (x>xmax) xmax=x; }
        avg += x;
        std += x*x;
        ct++;
    }

    avg /= ct;
    std = sqrt( (std/ct) - avg*avg );

    printf("Input file:     %s\tcolumn %d\n",argv[1],column);
    printf("# of data pts:  %d\n",ct);
    printf("Average:        %f\n",avg);
    printf("Std. Deviation: %f\n",std);
    printf("Min,Max values: %f,%f\n",xmin,xmax);
    exit(1);
}
