#include"mgmres.h"
#include"mgmres.c"
#include<stdio.h>
#include<math.h>

#define delta 0.1

/* do zapisu plikow*/
#include <sys/stat.h>
#include <unistd.h>

void create_directory(const char *path) {
    #if defined(_WIN32) || defined(_WIN64)
        _mkdir(path);
    #else
        mkdir(path, 0777);
    #endif
}
/*endif*/

int set_l(int i, int j, int nx){
    return i + j*(nx+1);
}//wzor (11)

int get_j(int nx, int l){
    return floor(l/(nx+1));
}//wzor (12)

int get_i(int nx, int l){
    return l - get_j(nx, l)*(nx+1);
}//wzor (13)

double rho_1(int i, int j, int nx, int ny){
    return 0.0;
}//w przykladach realizowanych w tym pliku rho_1 oraz rho_2 sa rowne 0

double rho_2(int i, int j, int nx, int ny){
    return 0.0;
}//w przykladach realizowanych w tym pliku rho_1 oraz rho_2 sa rowne 0

double eps_l(int eps1, int eps2, int l, int nx){
    if(get_i(nx, l) >= (nx/2)){
        return (double)eps1;
    }else{
        return (double)eps2;
    }
}//wzor (21)

int macierz_A(int N, int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4, double *a, int *ja, int *ia, double *b){
    int brzeg; //jesli 1 to jestesmy na brzegu, 0 to srodek
    double vb; //potencjal na brzegu
    int k = -1;

    if(nx == 4){
        create_directory("wypelnienie");
        chdir("wypelnienie");
        FILE* b_vec = fopen("wektor_b.txt", "w");
        fprintf(b_vec, "#l \t i_l \t j_l \t b[l]\n");
        fclose(b_vec);
    }//zapis do plikow wektora b przy nx=ny=4

    for(int l=0; l<N; l++){
        brzeg = 0;
        vb = 0;

        if(get_i(nx, l) == 0){ //lewy brzeg i = 0
            brzeg = 1;
            vb = V1;
        }
        if(get_i(nx, l) == nx){ //prawy brzeg i = nx
            brzeg = 1;
            vb = V3;
        }
        if(get_j(nx, l) == 0){ //dol
            brzeg = 1;
            vb = V4;
        }
        if(get_j(nx, l)==ny){ //gora
            brzeg = 1;
            vb = V2;
        }

        /*wypelnienie wektora wyrazow wolnych*/
        b[l] = -(rho_1(get_i(nx, l), get_j(nx, l), nx, ny)+rho_2(get_i(nx, l), get_j(nx, l), nx, ny));

        if(brzeg == 1){
            b[l] = vb;
        }

        if(nx==4){
            FILE* b_vect = fopen("wektor_b.txt", "a");
            fprintf(b_vect, "%d \t %d \t %d \t %lf \n", l, get_i(nx, l), get_j(nx, l), b[l]);
            fclose(b_vect);
        }//zapis do pliku

        /*wypelnienie elementow macierzy A*/
        ia[l] = -1; //wskaznik pierwszego elementu w wierszu

        if(l-nx-1 >= 0 && brzeg == 0){ //lewa skrajna diagonala
            k++;
            if(ia[l]<0){
                ia[l] = k;
            }
            a[k] = eps_l(eps1, eps1, l, nx)/(delta*delta); //a[l][l-nx-1], wzor (16)
            ja[k] = l - nx - 1;
        }
        if(l-1 >= 0 && brzeg == 0){ //poddiagonala
            k++;
            if(ia[l]<0){
                ia[l] = k;
            }
            a[k] = eps_l(eps1, eps2, l, nx)/(delta*delta); //a[l][l-1], wzor (17)
            ja[k] = l - 1;
        }
        //diagonala
        k++;
        if(ia[l]<0){
            ia[l] = k;
        }
        if(brzeg == 0){
            a[k] = -(2 * eps_l(eps1, eps2, l, nx) + eps_l(eps1, eps2, l+1, nx) + eps_l(eps1, eps2, l+nx+1, nx))/(delta*delta);//a[l][l], wzor (18)
        }else{
            a[k] = 1;
        }
        ja[k] = l;

        if(l<N && brzeg == 0){//naddiagonala
            k++;
            a[k] = eps_l(eps1, eps2, l+1, nx)/(delta*delta); //a[l][l+1], wzor (19)
            ja[k] = l+1;
        }
        if(1<(N-nx-1) && brzeg == 0){ //skrajna prawa diagonala
            k++;
            a[k] = eps_l(eps1, eps2, l+nx+1, nx)/(delta*delta);//a[l][l+nx+1], wzor (20)
            ja[k] = l + nx +1;
        }
    }
    int nz_num = k+1; //ilosc niezerowych elementow
    ia[N] = nz_num;

    if(nx==4){
        FILE *A = fopen("Macierz_A.txt", "w");
        fprintf(A, "#k \t a[k] \n");

        for(int k=0; k<5*N; k++){
            fprintf(A, "%d \t %lf \n", k, a[k]);
        }
        fclose(A);
        chdir("..");
    }//zapis do pliku macierzy A

    return nz_num;
}

void metoda_algebraiczna(int N, int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4, double *a, int *ja, int *ia, double *b, double *V){
    chdir("files");
    char filename[256];
    sprintf(filename, "%s%d%s", "map_nx_ny_", nx, ".txt");
    FILE* map = fopen(filename, "w");

    int nz_num = macierz_A(N, nx, ny, eps1, eps2, V1, V2, V3, V4, a, ja, ia, b); //wypelnienie macierzy

    /*parametry zgodnie z trescia zadania*/
    int itr_max = 500;
    int mr = 500;
    double tol_abs = pow(10, -8);
    double tol_rel = pow(10, -8);

    /*uzycie danej nam biblioteki*/
    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);

    for(int l=0; l<N; l++){
        int i = get_i(nx, l);
        int j = get_j(nx, l);
        fprintf(map, "%lf %lf %lf \n", delta*i,  delta*j, V[l]);
    }
    fclose(map);
    chdir("..");
}

int main(){
    create_directory("files");

    /*ustalenie parametrow poczatkowych*/
    int nx = 4;
    int ny = 4;
    int eps1 = 1;
    int eps2 = 1;
    int V1 = 10;
    int V2 = -10;
    int V3 = 10;
    int V4 = -10;

    int N = (nx+1) * (ny+1);

    /*definicja werktorow*/
    double *a = malloc(sizeof(double) * 5 * N);
    int *ja = malloc(sizeof(int) * 5 * N);
    int *ia = malloc(sizeof(int) * (N+1));
    for(int i=0; i<=N; i++){
        ia[i] = -1;
    }
    double *b = malloc(sizeof(double) * N);
    double *V = malloc(sizeof(double) * N);

    metoda_algebraiczna(N, nx, ny, eps1, eps2, V1, V2, V3, V4, a, ja, ia, b, V);

    /*podpunkt a*/
    nx = 50;
    ny = 50;
    N = (nx+1) * (ny+1);

    a = realloc(a, sizeof(double) * 5 * N);
    ja = realloc(ja, sizeof(int) * 5 * N);
    ia = realloc(ia, sizeof(int) * (N+1));
    b = realloc(b, sizeof(double)*N);
    V = realloc(V, sizeof(double)*N);

    metoda_algebraiczna(N, nx, ny, eps1, eps2, V1, V2, V3, V4, a, ja, ia, b, V);

    /*podpunkt b*/
    nx = 100;
    ny = 100;
    N = (nx+1) * (ny+1);

    a = realloc(a, sizeof(double) * 5 * N);
    ja = realloc(ja, sizeof(int) * 5 * N);
    ia = realloc(ia, sizeof(int) * (N+1));
    b = realloc(b, sizeof(double)*N);
    V = realloc(V, sizeof(double)*N);

    metoda_algebraiczna(N, nx, ny, eps1, eps2, V1, V2, V3, V4, a, ja, ia, b, V);

    /*podpunkt c*/
    nx = 200;
    ny = 200;
    N = (nx+1) * (ny+1);

    a = realloc(a, sizeof(double) * 5 * N);
    ja = realloc(ja, sizeof(int) * 5 * N);
    ia = realloc(ia, sizeof(int) * (N+1));
    b = realloc(b, sizeof(double)*N);
    V = realloc(V, sizeof(double)*N);

    metoda_algebraiczna(N, nx, ny, eps1, eps2, V1, V2, V3, V4, a, ja, ia, b, V);

    free(a);
    free(ja);
    free(ia);
    free(b);
    free(V);
}
