#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//definicja wartosci parametrow
#define epsilon 1.0
#define delta 0.1
#define nx 150
#define ny 100
#define V1 10.0
#define V2 0.0
#define xmax delta * nx
#define ymax delta * ny
#define sigma_x 0.1 * xmax
#define sigma_y 0.1 * ymax
#define omega_G_val_1 0.6
#define omega_G_val_2 1.0
#define TOL 1e-8

#include <sys/stat.h>
#include <unistd.h>

void create_directory(const char *path) {
    #if defined(_WIN32) || defined(_WIN64)
        _mkdir(path);
    #else
        mkdir(path, 0777);
    #endif
}

double rho_1(int x, int y) {
    return exp(-pow((x - 0.35 * xmax) / sigma_x, 2) - pow((y - 0.5 * ymax) / sigma_y, 2));
}//wzor (19)

double rho_2(int x, int y) {
    return -exp(-pow((x - 0.65 * xmax) / sigma_x, 2) - pow((y - 0.5 * ymax) / sigma_y, 2));
}//wzor (20)

double rho(int x, int y) {
    return rho_1(delta * x, delta * y) + rho_2(delta * x, delta * y);
}

void rho_tab(int x, int y, double r[x][y]){
    for(int i=0; i<x; i++){
        for(int j=0; j<y; j++){
            r[i][j] = rho(i, j);
        }
    }
}//calosc trzymana w tablicy dla latwiejszego dostepu danych

double S_calka(double V[nx+1][ny+1], double rho[nx+1][ny+1]) {
    double s = 0;

    for (int i = 0; i <= nx - 1; i++) {
        for (int j = 0; j <= ny - 1; j++) {
            double sum = delta*delta * (
                0.5 * pow((V[i + 1][j] - V[i][j]) / delta, 2) +
                0.5 * pow((V[i][j + 1] - V[i][j]) / delta, 2) -
                rho[i][j] * V[i][j]
            );
            s += sum;
        }
    }
    return s;
}//wyliczenie calki funkcjonalnej - wzor (16), wyliczone za pomoca wzoru (17)

void relaksacja_globalna(double omega_G, double rho_tab[nx+1][ny+1], const char* folder, const char* filename, const char* map_name, const char* err_name){
    printf("relaksacja globalna, omega = %f: \n", omega_G);

    char filepath[256];
    create_directory(folder);
    chdir(folder);

    double V_s[nx+1][ny+1] = {0.0}; //inicjalizacja tablic V_s
    double V_n[nx+1][ny+1] = {0.0}; //oraz V_n
    for (int i = 0; i <= nx; i++){
        V_s[i][0] = V1;
        V_s[i][ny] = V2;
        V_n[i][0] = V1;
        V_n[i][ny] = V2;
    }//warunki brzegowe

    int iteracja = 0; //licznik

    sprintf(filepath, "%s/%s", folder, filename);
    FILE* file = fopen(filename, "w");
    sprintf(filepath, "%s/%s", folder, map_name);
    FILE* map = fopen(map_name, "w");
    sprintf(filepath, "%s/%s", folder, err_name);
    FILE* err = fopen(err_name, "w");

    double S_it = S_calka(V_s, rho_tab); //inicjalizacja wartosci obecnej calki

    while (1) {
        iteracja++;

        for (int i = 1; i <= nx - 1; i++) {
            for (int j = 1; j <= ny - 1; j++) {
                V_n[i][j] = 0.25 * (
                    V_s[i + 1][j] + V_s[i - 1][j] +
                    V_s[i][j + 1] + V_s[i][j - 1] +
                    ((delta * delta) / epsilon) * rho_tab[i][j]); //wyznaczenie elementow nowej tablicy - wzor (9)
            }
        }
        for (int j = 1; j <= ny - 1; j++) {
            V_n[0][j] = V_n[1][j];
            V_n[nx][j] = V_n[nx - 1][j];
        }//WB Neumanna - wzor (10) i (11)

        for (int i = 0; i <= nx; i++) {
            for (int j = 0; j <= ny; j++) {
                V_s[i][j] = (1.0 - omega_G) * V_s[i][j] + omega_G * V_n[i][j]; //wymieszanie rozwiazan - wzor (12)
            }
        }

        double S_it_1 = S_it;
        S_it = S_calka(V_s, rho_tab); // nowe wyliczenie wartosci calki

        fprintf(file, "%d %.8lf\n", iteracja, S_it);

        if(fabs((S_it-S_it_1)/S_it_1)<TOL){
            break; //sprawdzenie warunku stopu - wzor (18)
        }

        if(iteracja%1000 == 1){
            printf("liczenie ... %d iteracja\n", iteracja);
        }
    }
    printf("zakonczono - %d iteracji\n", iteracja);
    
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            fprintf(map, "%.8lf %.8lf %.8lf \n", delta * i, delta * j, V_s[i][j]); //dane do mapy zrelaksowanego potencjalu
        }
    }

    for(int i=0; i<=nx; i++){
        for(int j=0; j<=nx; j++){
            double delta_smol = 0;
            if(i > 0 && i < nx && j > 0 && j < ny){
                delta_smol = (V_s[i+1][j] - 2 * V_s[i][j] + V_s[i-1][j])/(delta*delta) +
                (V_s[i][j+1] - 2 * V_s[i][j] + V_s[i][j-1])/(delta*delta) + (rho_tab[i][j]/epsilon); //dane do mapy bledu
            }
            fprintf(err, "%.8lf %.8lf %.8lf \n", delta * i, delta * j, delta_smol);
        }
    }

    fclose(file);
    fclose(map);
    fclose(err);
    chdir("..");
}

int main() {
    const char *folder = "files";

    double rho_all_values[nx+1][ny+1];
    rho_tab(nx+1, ny+1, rho_all_values);

    relaksacja_globalna(omega_G_val_1, rho_all_values, folder, "glob_1.txt", "map_omega_1.txt", "err_1.txt");
    relaksacja_globalna(omega_G_val_2, rho_all_values, folder, "glob_2.txt", "map_omega_2.txt", "err_2.txt");

    return 0;
}