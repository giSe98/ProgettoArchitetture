/**************************************************************************************
* 
* CdL Magistrale in Ingegneria Informatica
* Corso di Architetture e Programmazione dei Sistemi di Elaborazione - a.a. 2020/21
* 
* Progetto dell'algoritmo di Regressione
* in linguaggio assembly x86-32 + SSE
* 
* Fabrizio Angiulli, aprile 2019
* 
**************************************************************************************/

/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib32gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf64 regression64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o regression64.o regression46c.c -o regression64c -lm && ./regression64c $pars
* 
* oppure
* 
* ./runregression64
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

typedef struct {
    MATRIX x; //data set
    VECTOR y; //label set
    MATRIX xast; //data set convertito
    int n; //numero di punti del data set
    int d; //numero di dimensioni del data set    
    int k; //dimensione del batch
    int degree; //grado del polinomio
    type eta; //learning rate
    //STRUTTURE OUTPUT
    VECTOR theta; //vettore dei parametri
    int t; //numero di parametri, dimensione del vettore theta
    int iter; //numero di iterazioni
    int adagrad; //accelerazione adagrad
    int silent; //silenzioso
    int display; //stampa risultati
} params;

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (double*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (double**).
* 
* 	In entrambi i casi il candidato dovrà inoltre scegliere se memorizzare le
* 	matrici per righe (row-major order) o per colonne (column major-order).
*
* 	L'assunzione corrente è che le matrici siano in row-major order.
* 
*/

void* get_block(int size, int elements) { 
    return _mm_malloc(elements*size,32); 
}

void free_block(void* p) { 
    _mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
    return (MATRIX) get_block(sizeof(type),rows*cols);
}

void dealloc_matrix(MATRIX mat) {
    free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
* 
*****************************************************************************
*	Se lo si ritiene opportuno, è possibile cambiare la codifica in memoria
* 	della matrice. 
*****************************************************************************
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
    FILE* fp;
    int rows, cols, status, i;
    
    fp = fopen(filename, "rb");
    
    if (fp == NULL){
        printf("'%s': bad data file name!\n", filename);
        exit(0);
    }
    
    status = fread(&cols, sizeof(int), 1, fp);
    status = fread(&rows, sizeof(int), 1, fp);
    
    MATRIX data = alloc_matrix(rows,cols);
    status = fread(data, sizeof(type), rows*cols, fp);
    fclose(fp);
    
    *n = rows;
    *k = cols;
    
    return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di righe (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di colonne (M) --> numero intero a 32 bit
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
*/
void save_data(char* filename, void* X, int n, int k) {
    FILE* fp;
    int i;
    fp = fopen(filename, "wb");
    if(X != NULL){
        fwrite(&k, 4, 1, fp);
        fwrite(&n, 4, 1, fp);
        for (i = 0; i < n; i++) {
            fwrite(X, sizeof(type), k, fp);
            //printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
            X += sizeof(type)*k;
        }
    }
    fclose(fp);
}

// PROCEDURE ASSEMBLY

extern void prodottoScalare(double* a, double* b, int dim,double* pScalare);
extern void aggiornaTheta(double* theta,double* x, double fattore,int dim);
extern void aggiornaThetaAdagrad(double* theta, double* x, double fattore, int dim, MATRIX G, type pScalare);


// Returns value of Binomial Coefficient C(n, k) 
int binomialCoeff(int n, int k) { 
    int res = 1, i; 
  
    // Since C(n, k) = C(n, n-k) 
    if (k > n - k) 
        k = n - k; 
  
    // Calculate value of 
    // [n * (n-1) --- (n-k+1)] / [k * (k-1) ---- 1] 
    for (i = 0; i < k; ++i) { 
        res *= (n - i); 
        res /= (i + 1); 
    } 
  
    return res; 
}

void estendi(params* input,int offset_xast, int offset_x){
    //INIZIALIZZAZIONE PRIMI LENX+1 ELEMENTI
    input->xast[offset_xast]=1;
    for(int i = 0; i<input->d; ++i)
        input->xast[i+1+offset_xast]=input->x[i+offset_x];
    //START DA GRADO 2
    int* aux=malloc(input->degree *sizeof(int));
	int ind=input->d+1;
	int i=2;
	for (;i<=input->degree;++i) { 
		int temp=0;
		while(temp<i) {
			aux[temp]=0;
			temp++;
		}
		temp--;
		input->xast[ind+offset_xast]=pow(input->x[offset_x],i);
		ind++;
		while(temp>=1 && input->d != 1) {
			if(aux[temp]==input->d-1) {
				temp--;
				aux[temp]++;
				if(aux[temp]!=input->d-1)
					while(temp<i-1) {
						aux[temp+1]=aux[temp];
						temp++;
					}
			}
			else
				aux[temp]++;
				int k=0;
				input->xast[ind+offset_xast]=1;
				for(;k<i;k++) {
					input->xast[ind+offset_xast]*=input->x[aux[k]+offset_x];
				}
			ind++;
		}
	}
}

void convert_data(params* input){
    int upper = binomialCoeff(input->degree+input->d, input->d);
    input->t=upper;
    input->xast =alloc_matrix(input->n,upper);
    int i;
    //POSSIBILE LOOP-UNROLLING / CACHE BLOCKING -> LAVORANDO PER RIGHE INDIPENDENTI
    for(i = 0; i<input->n; ++i){
        //mettere direttiva openmp
        estendi(input,i*upper,i*input->d);
    }
}

void calcoloNuovoThetaAdagrad(int j, int v, MATRIX G, params* input){
    int i;
    type pScalare;
    type fattore = input->eta/v;
    for(i=j;i<j+v;++i){
        prodottoScalare(input->theta,&input->xast[i*input->t],input->t,&pScalare);
        aggiornaThetaAdagrad(input->theta,&input->xast[i*input->t],fattore,input->t,G,(pScalare- input->y[i]));
    }

}

void calcoloNuovoTheta(int j, int v, params* input){
    int i;
    type pScalare,fattore, param = input->eta/v;
    for(i=j; i<j+v; ++i){
        prodottoScalare(input->theta,&input->xast[i*input->t],input->t,&pScalare);
        fattore = (pScalare- input->y[i]) * param;
        aggiornaTheta(input->theta,&input->xast[i*input->t],fattore,input->t);
    }
}
void sgd(params* input){
   input->theta = alloc_matrix(1,input->t); //creazione vettore theta (parametri)
    int iter=0,i,j;
    MATRIX G;
    for(i=0;i<input->t;i++)
        input->theta[i]=0;
    if(input->adagrad){
        G = alloc_matrix(1,input->t);
	for(i=0;i<input->t;i++){
		G[i]=0;
	}
    }
    while(iter<input->iter){
        for(j=0;j<input->n;j+=input->k){
            int v = input->k;
            if(input->n-j<input->k)   v = input->n-j;
            //calcolo nuovo vettore theta
            
            if(input->adagrad)
                calcoloNuovoThetaAdagrad(j, v, G, input);
            else
                calcoloNuovoTheta(j, v, input);
        }
        iter+=1;
    }
}

int main(int argc, char** argv) {

    char fname[256];
    char* dsname;
    char* filename;
    int i, j, k;
    clock_t t;
    float time;
    int yd = 1;
    
    //
    // Imposta i valori di default dei parametri
    //

    params* input = malloc(sizeof(params));
    
    input->x = NULL;
    input->y = NULL;
    input->xast = NULL;
    input->n = 0;
    input->d = 0;
    input->k = -1;
    input->degree = -1;
    input->eta = -1;
    input->iter = -1;
    input->adagrad = 0;
    input->theta = NULL;
    input->t = 0;
    input->adagrad = 0;
    input->silent = 0;
    input->display = 0;

    //
    // Visualizza la sintassi del passaggio dei parametri da riga comandi
    //

    if(argc <= 1){
        printf("%s D -batch <k> -degree <deg> -eta <eta> -iter <it> [-adagrad]\n", argv[0]);
        printf("\nParameters:\n");
        printf("\tD: il nome del file, estensione .data per i dati x, estensione .labels per le etichette y\n");
        printf("\t-batch <k>: il numero di campini nel batch\n");
        printf("\t-degree <deg>: il grado del polinomio\n");
        printf("\t-eta <eta>: il learning rate\n");
        printf("\t-iter <it>: il numero di iterazioni\n");
        printf("\t-adagrad: l'acceleratore AdaGrad\n");
        exit(0);
    }
    
    //
    // Legge i valori dei parametri da riga comandi
    //
    
    int par = 1;
    while (par < argc) {
        if (par == 1) {
            filename = argv[par];
            par++;
        } else if (strcmp(argv[par],"-s") == 0) {
            input->silent = 1;
            par++;
        } else if (strcmp(argv[par],"-d") == 0) {
            input->display = 1;
            par++;
        } else if (strcmp(argv[par],"-batch") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing batch dimension value!\n");
                exit(1);
            }
            input->k = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-degree") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing degree value!\n");
                exit(1);
            }
            input->degree = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-eta") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing eta value!\n");
                exit(1);
            }
            input->eta = atof(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-iter") == 0) {
            par++;
            if (par >= argc) {
                printf("Missing iter value!\n");
                exit(1);
            }
            input->iter = atoi(argv[par]);
            par++;
        } else if (strcmp(argv[par],"-adagrad") == 0) {
            input->adagrad = 1;
            par++;
        } else{
            printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
            par++;
        }
    }
    
    //
    // Legge i dati e verifica la correttezza dei parametri
    //
    
    if(filename == NULL || strlen(filename) == 0){
        printf("Missing input file name!\n");
        exit(1);
    }

    dsname = basename(strdup(filename));
    sprintf(fname, "%s.data", filename);
    input->x = load_data(fname, &input->n, &input->d);
    sprintf(fname, "%s.labels", filename);
    input->y = load_data(fname, &input->n, &yd);

    if(input->k < 0){
        printf("Invalid value of batch dimension parameter!\n");
        exit(1);
    }
    
    if(input->degree < 0){
        printf("Invalid value of degree parameter!\n");
        exit(1);
    }
    
    if(input->eta < 0){
        printf("Invalid value of eta parameter!\n");
        exit(1);
    }
    
    if(input->iter < 0){
        printf("Invalid value of iter parameter!\n");
        exit(1);
    }
    
    //
    // Visualizza il valore dei parametri
    //
    
    if(!input->silent){
        printf("Input data name: '%s.data'\n", filename);
        printf("Input label name: '%s.labels'\n", filename);
        printf("Data set size [n]: %d\n", input->n);
        printf("Number of dimensions [d]: %d\n", input->d);
        printf("Batch dimension: %d\n", input->k);
        printf("Degree: %d\n", input->degree);
        printf("Eta: %f\n", input->eta);
        if(input->adagrad)
            printf("Adagrad enabled\n");
        else
            printf("Adagrad disabled\n");
    }
    
    //
    //prova(input);
    //

    //
    // Conversione Dati
    //
    
    t = clock();
    convert_data(input);
    t = clock() - t;
    time = ((float)t)/CLOCKS_PER_SEC;
    sprintf(fname, "%s.xast", dsname);
       
    if(!input->silent)
        printf("Conversion time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);
    
    //
    // Regressione
    //
    
    t = clock();
    sgd(input);
    t = clock() - t;
    time = ((float)t)/CLOCKS_PER_SEC;
    
    if(!input->silent)
        printf("Regression time = %.3f secs\n", time);
    else
        printf("%.3f\n", time);

    //
    // Salva il risultato di theta
    //
    
    if(!input->adagrad)
	    sprintf(fname, "%s.theta.sgd", dsname);
    else
	    sprintf(fname, "%s.theta.adagrad", dsname);
    save_data(fname, input->theta, input->t, 1);
    if(input->display){
        printf("theta: [");
        for(i=0; i<input->t-1; i++)
            printf("%f,", input->theta[i]);
        printf("%f]\n", input->theta[i]);
    }
    
    if(!input->silent)
        printf("\nDone.\n");

    return 0;
}
