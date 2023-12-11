#include<stdio.h>
#include<math.h>

#define delta 0.01
#define rho 1
#define mi 1
#define nx 200
#define ny 90
#define i_1 50
#define j_1 55
#define IT_MAX 20000

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


double y_j(double j){
    return delta * j;
}

double x_i(double i){
    return delta * i;
}

void brzeg_psi(double psi[nx+1][ny+1], double Q_we){
    double y_j1 = y_j(j_1);
    double y_ny = y_j(ny);

    double Q_wy = Q_we * (pow(y_ny, 3) - pow(y_j1, 3) - 3 * y_j1 * pow(y_ny, 2) + 3 * pow(y_j1, 2) * y_ny)/(pow(y_ny, 3));
    double y;

    //brzeg A (wejscie)
    for(int j = j_1; j<=ny; j++){
        y = delta * j;
        psi[0][j] = (Q_we/(2*mi)) * ((pow(y, 3)/3) - (pow(y, 2)/2) * (y_j1 + y_ny) + y * y_j1 * y_ny);
    }

    //brzeg C (wyjscie)
    for(int j=0; j<=ny; j++){
        y = delta * j;
        psi[nx][j] = (Q_wy/(2*mi)) * ((pow(y, 3)/3) - (pow(y, 2)/2)*y_ny) + (Q_we*pow(y_j1, 2)*(-y_j1+3*y_ny))/(12*mi);
    }

    //brzeg B
    for(int i=1; i<=nx-1; i++){
        psi[i][ny] = psi[0][ny];
    }

    //brzeg D
    for(int i=i_1; i<=nx-1; i++){
        psi[i][0] = psi[0][j_1];
    }

    //brzeg E
    for(int j=1; j<=j_1; j++){
        psi[i_1][j] = psi[0][j_1];
    }

    //brzeg F
    for(int i=1; i<=i_1; i++){
        psi[i][j_1] = psi[0][j_1];
    }
}

void brzeg_dzeta(double dzeta[nx+1][ny+1], double psi[nx+1][ny+1], double Q_we){
    double y_j1 = y_j(j_1);
    double y_ny = y_j(ny);

    double Q_wy = Q_we * (pow(y_ny, 3) - pow(y_j1, 3) - 3 * y_j1 * pow(y_ny, 2) + 3 * pow(y_j1, 2) * y_ny)/(pow(y_ny, 3));
    double y;

    //brzeg A (wejscie)
    for(int j=j_1; j<=ny; j++){
        y = delta*j;
        dzeta[0][j] = (Q_we/(2*mi)) * (2*y - y_j1 - y_ny);
    }

    //brzeg C (wyjscie)
    for(int j=0; j<=ny; j++){
        y = delta*j;
        dzeta[nx][j] = (Q_wy/(2*mi)) * (2*y - y_ny);
    }

    //brzeg B
    for(int i=1; i<=nx-1; i++){
        dzeta[i][ny] = (2.0/(delta*delta)) * (psi[i][ny-1] - psi[i][ny]);
    }

    //brzeg D
    for(int i=i_1+1; i<=nx-1; i++){
        dzeta[i][0] = (2.0/(delta*delta)) * (psi[i][1] - psi[i][0]);
    }

    //brzeg E
    for(int j=1; j<=j_1-1; j++){
        dzeta[i_1][j] = (2.0/(delta*delta)) * (psi[i_1+1][j] - psi[i_1][j]);
    }

    //brzeg F
    for(int i=1; i<=i_1; i++){
        dzeta[i][j_1] = (2.0/(delta*delta)) * (psi[i][j_1+1] - psi[i][j_1]);
    }

    //wierzcholek E/F
    dzeta[i_1][j_1] = 0.5 * (dzeta[i_1-1][j_1]+dzeta[i_1][j_1-1]);

}

void relaksacja_NS(double dzeta[nx+1][ny+1], double psi[nx+1][ny+1], double Q_we){
    create_directory("files");
    chdir("files");

    char filename[256];
    sprintf(filename, "Q_%d.txt", (int)Q_we);
    FILE * f = fopen(filename, "w");

    double u[nx+1][ny+1];
    double v[nx+1][ny+1];

    int omega;
    double Gamma;

    brzeg_psi(psi, Q_we);

    for(int it=1; it<=IT_MAX; it++){
        if(it<2000){
            omega = 0;
        }else{
            omega = 1;
        }

        for(int i=1; i<=nx-1; i++){
            for(int j=1; j<=ny-1; j++){
                if((i>i_1) || (j>j_1)){
                    psi[i][j] = 0.25 * (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] + psi[i][j-1] - delta*delta*dzeta[i][j]);
                    dzeta[i][j] = 0.25 * (dzeta[i+1][j] + dzeta[i-1][j] + dzeta[i][j+1] + dzeta[i][j-1]);
                    double var = (psi[i][j+1] - psi[i][j-1])*(dzeta[i+1][j] - dzeta[i-1][j]) - (psi[i+1][j] - psi[i-1][j])*(dzeta[i][j+1] - dzeta[i][j-1]);
                    dzeta[i][j] = dzeta[i][j] - omega*(rho/(16*mi))*var;
                }
                u[i][j] = (psi[i][j+1] - psi[i][j-1])/(2*delta);
                v[i][j] = -(psi[i+1][j] - psi[i-1][j])/(2*delta);
            }
        }

        brzeg_dzeta(dzeta, psi, Q_we);

        Gamma = 0.0;
        for(int i=1; i<=nx-1; i++){
            int j_2 = j_1+2;
            double tmp = psi[i+1][j_2] + psi[i-1][j_2] + psi[i][j_2+1] + psi[i][j_2-1] - 4*psi[i][j_2] - delta*delta*dzeta[i][j_2];
            Gamma = Gamma + tmp;
        }
        printf("IT = %d, gamma = %.10lf \n", it, Gamma);
    }

    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            double x = delta*i;
            double y = delta*j;
            fprintf(f, "%lf %lf %g %g %g %g\n", x, y, psi[i][j], dzeta[i][j], u[i][j], v[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
    chdir("..");
}

int main(){
    double psi[nx+1][ny+1];
    double dzeta[nx+1][ny+1];
    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            psi[i][j] = 0.0;
            dzeta[i][j] = 0.0;
        }
    }
    double Q = -1000;

    relaksacja_NS(dzeta, psi, Q);

    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            psi[i][j] = 0.0;
            dzeta[i][j] = 0.0;
        }
    }
    Q = -4000;

    relaksacja_NS(dzeta, psi, Q);

    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            psi[i][j] = 0.0;
            dzeta[i][j] = 0.0;
        }
    }
    Q = 4000;

    relaksacja_NS(dzeta, psi, Q);

    create_directory("img");

    return 0;
}