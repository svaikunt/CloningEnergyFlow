typedef struct{
    char *nam;
    FILE *ptr;
} FileWithName;

#define COMMENTS_LEN 10000
#define LINE_LEN 200

int WhichBin_Float(float x, float *b, int n);
int WhichBin_Double(double x, double *b, int n);
FILE *Open(char *filename, char *mode);
char *ReadComments(FILE *fp_in);
char *ReadAndWriteComments(FILE *fp_in, FILE *fp_out);
void Printf(char *string);
void ShowTime();
void PrintCommandLine(FILE *f, int argc, char **argv);
int **AllocateRectangularMatrix_Int(int dim1, int dim2);
float **AllocateRectangularMatrix_Float(int dim1, int dim2);
double **AllocateRectangularMatrix_Double(int dim1, int dim2);
char **AllocateRectangularMatrix_Char(int dim1, int dim2);
int CountLinesRemainingInFile(FILE *fp);
int CountLengthToEndOfCurrentLine(FILE *fp);
int Amino(char s);
char AminoChar(int i);
FILE *DrawLine(FILE *fp);
FILE *TimeStamp(FILE *fp);
int IsLetter(char c);
char UpperCase(char c);
char LowerCase(char c);
char *UpperCase_String(char *s);
char *LowerCase_String(char *s);
int RemoveBlanks(char *s);

int amino_util[]={1,0,2,3,4,5,6,7,8,0,9,10,11,12,0,13,14,15,16,17,0,18,19,0,20,0};

/*----------------------------------------------------------------*/

int WhichBin_Float(float x, float *b, int n)

{
    int i,low=0,high=n;

    if (x > b[n] || x < b[0]) return -1;

    while (high-low > 1){
        i = (high+low)/2;
        if (x <= b[i]) high = i;
        else low = i;
    }
    return low;
}

/*----------------------------------------------------------------*/

int WhichBin_Double(double x, double *b, int n)

{
    int i,low=0,high=n;

    if (x > b[n] || x < b[0]) return -1;

    while (high-low > 1){
        i = (high+low)/2;
        if (x <= b[i]) high = i;
        else low = i;
    }
    return low;
}

/*----------------------------------------------------------------*/

FILE *Open(char *filename, char *mode)

{
    FILE *fp;

    if ((fp=fopen(filename,mode))==NULL){
        printf("Can't open %s !\n",filename);
        exit(1);
    }

    return fp;
}

/*----------------------------------------------------------------*/

char *ReadComments(FILE *fp_in)

{
    long current;
    char line[LINE_LEN],*comments;

    comments = malloc(COMMENTS_LEN*sizeof(char));
    *comments='\0';

    current = ftell(fp_in);
    fgets(line,LINE_LEN,fp_in);
    while (*line=='%' || *line=='#'){
        strcat(comments,line);
        current = ftell(fp_in);
        fgets(line,LINE_LEN,fp_in);
    }
    fseek(fp_in,current,SEEK_SET);

    return comments;
}

/*----------------------------------------------------------------*/

char *ReadAndWriteComments(FILE *fp_in, FILE *fp_out)

{
    long current;
    char line[LINE_LEN],*comments;

    comments = malloc(COMMENTS_LEN*sizeof(char));
    *comments='\0';

    current = ftell(fp_in);
    fgets(line,LINE_LEN,fp_in);
    while (*line=='%' || *line=='#'){
        strcat(comments,line);
        current = ftell(fp_in);
        fgets(line,LINE_LEN,fp_in);
    }
    fseek(fp_in,current,SEEK_SET);

    fprintf(fp_out,"%s",comments);

    return comments;
}

/*----------------------------------------------------------------*/

void Printf(char *string)

{
    printf("%s",string);
    fflush(stdout);
}

/*----------------------------------------------------------------*/

void ShowTime()

{
    system("date +\"\t\t\t\t\t%a %b %d %Y %T\"");
    fflush(stdout);
}

/*----------------------------------------------------------------*/

void PrintCommandLine(FILE *f, int argc, char **argv)

{
    int i;

    fprintf(f,"#");
    for (i=0;i<argc;i++) fprintf(f," %s",argv[i]);
    fprintf(f,"\n");
    fflush(f);
}

/*----------------------------------------------------------------*/

int **AllocateRectangularMatrix_Int(int dim1, int dim2)

{
    int s,t;
    int **matrix;

    matrix = malloc(dim1*sizeof(int *));
    for (s=0;s<dim1;s++){
        matrix[s] = malloc(dim2*sizeof(int));
        for (t=0;t<dim2;t++) matrix[s][t] = 0;
    }

    return matrix;
}

/*----------------------------------------------------------------*/

float **AllocateRectangularMatrix_Float(int dim1, int dim2)

{
    int s,t;
    float **matrix;

    matrix = malloc(dim1*sizeof(float *));
    for (s=0;s<dim1;s++){
        matrix[s] = malloc(dim2*sizeof(float));
        for (t=0;t<dim2;t++) matrix[s][t] = 0.;
    }

    return matrix;
}

/*----------------------------------------------------------------*/

double **AllocateRectangularMatrix_Double(int dim1, int dim2)

{
    int s,t;
    double **matrix;

    matrix = malloc(dim1*sizeof(double *));
    for (s=0;s<dim1;s++){
        matrix[s] = malloc(dim2*sizeof(double));
        for (t=0;t<dim2;t++) matrix[s][t] = 0.;
    }

    return matrix;
}

/*----------------------------------------------------------------*/

char **AllocateRectangularMatrix_Char(int dim1, int dim2)

{
    int s,t;
    char **matrix;

    matrix = malloc(dim1*sizeof(char *));
    for (s=0;s<dim1;s++){
        matrix[s] = malloc(dim2*sizeof(char));
        for (t=0;t<dim2;t++) matrix[s][t] = '\0';
    }

    return matrix;
}

/*----------------------------------------------------------------*/

int CountLinesRemainingInFile(FILE *fp)

{
    int nl=0;
    char c;
    long offset;

    offset = ftell(fp);
    while ( (c=fgetc(fp)) != EOF) if (c=='\n') nl++;
    fseek(fp,offset,SEEK_SET);
    return nl;
}


/*----------------------------------------------------------------*/

int CountLengthToEndOfCurrentLine(FILE *fp)

{
    int len=0;
    char c;
    long offset;

    offset = ftell(fp);
    c = fgetc(fp);
    while ( c!=EOF && c!='\n'){
        len++;
        c = fgetc(fp);
    }
    fseek(fp,offset,SEEK_SET);
    return len;
}

/*----------------------------------------------------------------*/

int Amino(char s)

{
    int i;

    i = s-'A';
    if (i>=0 && i<26) return amino_util[i];
    else return 0;
}

/*----------------------------------------------------------------*/

char AminoChar(int i)

{
    int j;
    if (i<1 || i>20) return 'X';
    for (j=0; j<26; j++) if (amino_util[j]==i) return 'A'+j;
}

/*----------------------------------------------------------------*/

FILE *DrawLine(FILE *fp)

{
    int i,length=70;
    char *line,c='-';

    line = malloc((length+2)*sizeof(char));

    for (i=0;i<length;i++) line[i]=c;
    line[i++] = '\n';
    line[i++] = '\0';

    fprintf(fp,"%s",line);

    free(line);
    return fp;
}

/*----------------------------------------------------------------*/

FILE *TimeStamp(FILE *fp)

{
    time_t *tp;

    tp = malloc(sizeof(time_t));
    time(tp);
    fprintf(fp,"\t\t\t\t\t%s",asctime(localtime(tp)));
    fflush(fp);
    free(tp);
    return fp;
}

/*----------------------------------------------------------------*/

int IsLetter(char c)

{
    int u,l;

    u = c-'A';
    if (u>=0 && u<26) return u+1;

    l = c-'a';
    if (l>=0 && l<26) return -l-1;

    return 0;
}

/*----------------------------------------------------------------*/

char UpperCase(char c)

{
    if (c-'A'>=0 && c-'A'<26) return c;
    else if (c-'a'>=0 && c-'a'<26) return c+'A'-'a';
    else return c;
}

/*----------------------------------------------------------------*/

char LowerCase(char c)

{
    if (c-'A'>=0 && c-'A'<26) return c+'a'-'A';
    else if (c-'a'>=0 && c-'a'<26) return c;
    else return c;
}

/*----------------------------------------------------------------*/

char *UpperCase_String(char *s)

{
    int i=-1;

    while (s[++i]!='\0') s[i] = UpperCase(s[i]);
    return s;
}

/*----------------------------------------------------------------*/

char *LowerCase_String(char *s)

{
    int i=-1;

    while (s[++i]!='\0') s[i] = LowerCase(s[i]);
    return s;
}

/*----------------------------------------------------------------*/

int RemoveBlanks(char *s)

{
    int n_blanks=0,i=0,j;

    while (i<strlen(s)){
        if (s[i]==' '){
           for (j=i;j<strlen(s)-1;j++) s[j]=s[j+1];
           s[j]='\0';
           n_blanks++;
           }
        else i++;
    }
    return n_blanks;
}
