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
/*end*/

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
    double x_max = delta * nx;
    double y_max = delta * ny;
    double sigma = x_max/10.0;
    double x = i*delta;
    double y = j*delta;

    double var1 = pow((x - 0.25*x_max),2)/pow(sigma, 2);
    double var2 = pow((y - 0.5*y_max), 2)/pow(sigma, 2);

    return exp(-var1-var2);
}//wzor (25)

double rho_2(int i, int j, int nx, int ny){
    double x_max = delta * nx;
    double y_max = delta * ny;
    double sigma = x_max/10.0;
    double x = i*delta;
    double y = j*delta;

    double var1 = pow((x - 0.75*x_max),2)/pow(sigma, 2);
    double var2 = pow((y - 0.5*y_max), 2)/pow(sigma, 2);

    return (-1)*exp(-var1-var2);
}//wzor (26)

double eps_l(int eps1, int eps2, int l, int nx){
    double i = (double)get_i(nx, l);
    double var = (double)nx;
    var = var/2.0;
    if(i <= var){
        return (double)eps1;
    }else{
        return (double)eps2;
    }
}//wzor (21)

int macierz_A(int N, int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4, double *a, int *ja, int *ia, double *b){
    int brzeg; //jesli 1 to jestesmy na brzegu, 0 to srodek
    double vb; //potencjal na brzegu
    int k = -1;

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
        b[l] = -(rho_1(get_i(nx, l), get_j(nx, l), nx, ny) +rho_2(get_i(nx, l), get_j(nx, l), nx, ny));

        if(brzeg == 1){
            b[l] = vb;
        }

        /*wypelnienie elementow macierzy A*/
        ia[l] = -1; //wskaznik pierwszego elementu w wierszu
        double eps = 0.0;

        if(l-nx-1 >= 0 && brzeg == 0){ //lewa skrajna diagonala
            k++;
            if(ia[l]<0){
                ia[l] = k;
            }
            eps = eps_l(eps1, eps2, l, nx);
            a[k] = eps/(delta*delta); //a[l][l-nx-1]
            ja[k] = l - nx - 1;
        }
        if(l-1 >= 0 && brzeg == 0){ //poddiagonala
            k++;
            if(ia[l]<0){
                ia[l] = k;
            }
            eps = eps_l(eps1, eps2, l, nx);
            a[k] = eps/(delta*delta); //a[l][l-1]
            ja[k] = l - 1;
        }
        //diagonala
        k++;
        if(ia[l]<0){
            ia[l] = k;
        }
        if(brzeg == 0){
            eps = eps_l(eps1, eps2, l, nx);
            double eps_l1 = eps_l(eps1, eps2, l+1, nx);
            double eps_lnx1 = eps_l(eps1, eps2, l+nx+1, nx);
            a[k] = -(2 * eps + eps_l1 + eps_lnx1)/(delta*delta);//a[l][l]
        }else{
            a[k] = 1;
        }
        ja[k] = l;

        if(l<N && brzeg == 0){//naddiagonala
            k++;
            eps = eps_l(eps1, eps2, l+1, nx);
            a[k] = eps/(delta*delta); //a[l][l+1]
            ja[k] = l+1;
        }
        if(1<(N-nx-1) && brzeg == 0){ //skrajna prawa diagonala
            k++;
            eps = eps_l(eps1, eps2, l+nx+1, nx);
            a[k] = eps/(delta*delta);//a[l][l+nx+1]
            ja[k] = l + nx +1;
        }
    }
    int nz_num = k+1; //ilosc niezerowych elementow
    ia[N] = nz_num;

    return nz_num;
}

void metoda_algebraiczna(int N, int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4, double *a, int *ja, int *ia, double *b, double *V){
    chdir("files");
    char filename[256];
    sprintf(filename, "%s%d%s%d%s", "map_eps_", eps1, "_eps_", eps2, ".txt");
    FILE* map = fopen(filename, "w");
    
    int nz_num = macierz_A(N, nx, ny, eps1, eps2, V1, V2, V3, V4, a, ja, ia, b); //wypelnienie macierzy

    /*ustalenie wartosci parametrow zgodnie z trescia*/
    int itr_max = 500;
    int mr = 500;
    double tol_abs = pow(10, -8);
    double tol_rel = pow(10, -8);

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);

    for(int l=0; l<N; l++){
        int i = get_i(nx, l);
        int j = get_j(nx, l);
        fprintf(map, "%lf %lf %lf \n", delta*i, delta*j, V[l]);
    }//zapis mapy do pliku

    fclose(map);
    chdir("..");
}

int main(){
    create_directory("files");

    /*ustalenie parametrow poczatkowych*/
    int nx = 100;
    int ny = 100;
    int V1 = 0;
    int V2 = 0;
    int V3 = 0;
    int V4 = 0;

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
    for(int i=0; i<N; i++){
        V[i] = 0.0;
    }

    /*podpunkt a*/
    int eps1 = 1;
    int eps2 = 1;

    metoda_algebraiczna(N, nx, ny, eps1, eps2, V1, V2, V3, V4, a, ja, ia, b, V);

    /*podpunkt b*/
    for(int i=0; i<N; i++){
        V[i] = 0.0;
    }

    eps1 = 1;
    eps2 = 2;

    metoda_algebraiczna(N, nx, ny, eps1, eps2, V1, V2, V3, V4, a, ja, ia, b, V);

    /*podpunkt c*/
    for(int i=0; i<N; i++){
        V[i] = 0.0;
    }

    eps1 = 1;
    eps2 = 10;

    metoda_algebraiczna(N, nx, ny, eps1, eps2, V1, V2, V3, V4, a, ja, ia, b, V);

    free(a);
    free(ja);
    free(ia);
    free(b);
    free(V);
}
