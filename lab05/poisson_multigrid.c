#include<stdio.h>
#include<math.h>

//definiowanie wartosci parametrow dla zadania
#define delta 0.2
#define Nx 128
#define Ny 128
#define xMax (delta)*(Nx)
#define yMax (delta)*(Ny)
#define TOL 10e-8

//warunki brzegowe Dirichleta
#define VB1(y) (sin(M_PI * (y)/(yMax))) //wzor (17)
#define VB2(x) (-sin(2*M_PI * (x)/(xMax) )) //wzor (18)
#define VB3(y) (sin(M_PI * (y)/(yMax) )) //wzor(19)
#define VB4(x) (sin(2*M_PI * (x)/(xMax) )) //wzor(20)

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

double S_calka(double V[Nx + 1][Ny + 1], int k, double delta_val) {
    double s;
    for (int i = 0; i <= Nx - k; i += k) {
        for (int j = 0; j <= Ny - k; j += k) {
            double sum = (k * k * delta_val * delta_val / 2) * (
                pow((V[i+k][j] - V[i][j] + V[i+k][j+k] - V[i][j+k])/(2*k*delta_val), 2) +
                pow((V[i][j+k] - V[i][j] + V[i+k][j+k] - V[i+k][j])/(2*k*delta_val), 2)
            );
            s = s + sum;
        }
    }
    return s;
}//warunek stopu, calka wzor (14) -> suma wzor (15)

void relaksacja_wielosiatkowa(int k1){
    int k = k1;
    chdir("files");
    double V[Nx+1][Ny+1];

    for(int i=0; i<=Nx; i++){
        for(int j=0; j<=Ny; j++){
            V[i][j] = 0.0;
        }
    } //inicjalizacja siatki obliczeniowej

    //nakladanie warunkow brzegowych
    for(int j=0; j<=Ny; j++){
        double t = (double)j;
        V[0][j] = VB1(t*delta); //lewo
        V[Nx][j] = VB3(t*delta); //prawo
    }

    for(int i=0; i<=Nx; i++){
        double t = (double)i;
        V[i][Ny] = VB2(t*delta); //gora
        V[i][0] = VB4(t*delta); //dol
    }

    double S_it_1 = S_calka(V, k, delta); //poprzednia iteracja
    double S_it = S_it_1; //terazniejsza iteracja

    int iter_total = 0; //licznik iteracji dla calosci

    FILE* file_S = fopen("s_it.txt", "w");
    FILE* file_it = fopen("it.txt", "w");
    
    while(k>0){
        int iteracja = 0; //licznik iteracji dla jednego k

        while(1){
            iteracja++;
            iter_total++;

            for(int i=k; i<=Nx-k; i+=k){
                for(int j=k; j<=Ny-k; j+=k){
                    V[i][j] = 0.25 * (V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]); //wzor (7)
                }
            }

            S_it_1 = S_it;
            S_it = S_calka(V, k, delta); //obliczanie nowej wartosci calki

            fprintf(file_S, "%d %.8f\n", iter_total, S_it);

            double war = fabs((S_it - S_it_1)/S_it_1); //wyliczenie wartosci do porownania z TOL

            if(war < TOL){
                break;
            }//sprawdzenie warunku stopu
        }
        fprintf(file_it, "%d\n", iteracja);

        char filename[256];
        sprintf(filename, "%s%d%s", "map_k_", k,".txt");

        FILE* map = fopen(filename, "w");

        for(int i=0; i<=Nx; i+=k){
            for(int j=0; j<=Ny; j+=k){
                fprintf(map, "%.8lf %.8lf %.8lf \n", delta * i, delta * j, V[i][j]); //dane do mapy zrelaksowanego potencjalu
            }
        }

        fclose(map);


        ///
        /* zageszczenie siatki */
        ///

        for(int i=0; i<=Nx-k; i+=k){
            for(int j=0; j<=Ny-k; j+=k){
                V[i+k/2][j+k/2] = 0.25*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]); //wzor (8)
                if(i != Nx - k){ //uwaga 2 - warunki brzegowe wyznaczamy tylko raz
                    V[i+k][j+k/2] = 0.5 * (V[i+k][j] + V[i+k][j+k]); //wzor (9)
                }
                if(j != Ny - k){
                    V[i+k/2][j+k] = 0.5 * (V[i][j+k] + V[i+k][j+k]); //wzor (10)
                }
                if(j != 0){
                    V[i+k/2][j] = 0.5 * (V[i][j] + V[i+k][j]); //wzor (11)
                }
                if(i != 0){
                    V[i][j+k/2] = 0.5 * (V[i][j] + V[i][j+k]); //wzor (12)
                }
            }
        }

        k = k/2; //przejscie do nastepnej wartosci k
    }
    fclose(file_S);
    chdir("..");
}

int main(){
    create_directory("files");

    relaksacja_wielosiatkowa(16);
}