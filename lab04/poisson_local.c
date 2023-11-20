#include <stdio.h>
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
}//calosc trzymam od razu w tablicy dla latwiejszego dostepu do danych

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

void relaksacja_lokalna(double omega_L, double rho_tab[nx+1][ny+1], const char* folder, const char* filename){
    printf("relaksacja lokalna, omega = %f: \n", omega_L);

    char filepath[256];
    create_directory(folder);
    chdir(folder);

    double V[nx+1][ny+1] = {0.0};//inicjalizaca tablicy V
    for(int i=0; i<=nx; i++){
        V[i][0] = V1;
        V[i][ny] = V2;
    }//na calosci sa 0 oprocz warunkow brzegowych

    int iteracja = 0; //licznik
    double S_it = S_calka(V, rho_tab); //inicjalizacja wartosci calki

    sprintf(filepath, "%s/%s", folder, filename);
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Błąd otwarcia pliku %s\n", filename);
    }

    while (1) {
        iteracja++;

        for (int i = 1; i <= nx-1; i++) {
            for (int j = 1; j <= ny-1; j++) {
                V[i][j] = (1.0-omega_L) * V[i][j] + (omega_L/4.0)*(
                    V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] + (((delta*delta)/epsilon) * rho_tab[i][j])
                ); //wzor (13) - nowe wartosci tablicy V
            }
        }

        for (int j = 1; j <= ny-1; j++) {
            V[0][j] = V[1][j];
            V[nx][j] = V[nx - 1][j];
        } //WB von Neumanna - wzory (15) i (16)

        double S_it_1 = S_it;
        S_it = S_calka(V, rho_tab); //nowa wartosc calki

        fprintf(file, "%d %.8lf\n", iteracja, S_it);

        if(fabs((S_it-S_it_1)/S_it_1)<TOL){
            break; //sprawdzenie warunku stopu - wzor (18)
        }

        if(iteracja%500 == 1){
            printf("liczenie ... %d iteracja\n", iteracja);
        }
    }
    printf("zakonczono - %d iteracji\n", iteracja);
    fclose(file);
    chdir("..");
}

int main() {
    const char *folder = "files";

    double omega_L_values[4] = {1.0, 1.4, 1.8, 1.9}; //inicjalizaca zmiennych wartosci omega_L

    double rho_all_values[nx+1][ny+1];
    rho_tab(nx+1, ny+1, rho_all_values);

    relaksacja_lokalna(omega_L_values[0], rho_all_values, folder, "lok_1.txt");
    relaksacja_lokalna(omega_L_values[1], rho_all_values, folder, "lok_2.txt");
    relaksacja_lokalna(omega_L_values[2], rho_all_values, folder, "lok_3.txt");
    relaksacja_lokalna(omega_L_values[3], rho_all_values, folder, "lok_4.txt");

    return 0;
}