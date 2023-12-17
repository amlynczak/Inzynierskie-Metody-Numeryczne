#include<stdio.h>
#include<math.h>
#include <sys/stat.h>

/* wyznaczenie parametrow do obliczen*/
#define nx 400
#define ny 90
#define i1 200
#define i2 210
#define j1 50
#define delta 0.01
#define sigma (10 * delta)
#define xA 0.45
#define yA 0.45

#define IT_MAX 10000 
//ustalone tak, aby bylo widac trzy maksima na wykresie x_sr

void wczytaj_strumien(double psi[nx+1][ny+1]){ //wczytanie danych strumienia z pliku psi.dat
    FILE *f;
    f = fopen("psi.dat", "r");

    if(f==NULL){
        printf("nie otwarto pliku psi.dat\n");
    }

    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            int r_i, r_j;
            double r_psi;

            fscanf(f, "%d %d %lf", &r_i, &r_j, &r_psi);

            if(r_i == i && r_j == j){
                psi[i][j] = r_psi;
            }
        }
    }
    fclose(f);
}

void setVxy(double Vx[nx+1][ny+1], double Vy[nx+1][ny+1], double psi[nx+1][ny+1]){ //wyznaczanie pola predkosci (Vx oraz Vy)  - wzory (11) - (16)
    for(int i=1; i<=nx-1; i++){
        for(int j=1; j<=ny-1; j++){
            Vx[i][j] = (psi[i][j+1] - psi[i][j-1])/(2*delta);
            Vy[i][j] = -(psi[i+1][j] - psi[i-1][j])/(2*delta);
        }
    }//wzor (11) i (12)

    for(int i = i1; i<=i2; i++){
        for(int j = 0; j<=j1; j++){
            Vx[i][j] = 0;
            Vy[i][j] = 0;
        }
    }//wzor (13)

    for(int i=1; i<=nx-1; i++){
        Vx[i][0] = 0;
        Vy[i][ny] = 0;
    }//wzor (14)

    for(int j = 0; j<=ny; j++){
        Vx[0][j] = Vx[1][j];
        Vx[nx][j] = Vx[nx-1][j];
    }//wzory (15) i (16)

    FILE *fvx, *fvy;
    fvx = fopen("files/Vx.txt", "w");
    fvy = fopen("files/Vy.txt", "w");

    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            fprintf(fvx, "%lf %lf %lf \n", i*delta, j*delta, Vx[i][j]);
            fprintf(fvy, "%lf %lf %lf \n", i*delta, j*delta, Vy[i][j]);
        }
    }

    fclose(fvx);
    fclose(fvy);
}

double find_dt(double Vx[nx+1][ny+1], double Vy[nx+1][ny+1]){
    double v_max = 0;
    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            double v_mod = sqrt(pow(Vx[i][j], 2)+pow(Vy[i][j], 2));
            if(v_mod > v_max) v_max = v_mod;
        }
    }
    double dt = delta/(4*v_max);
    return dt;
}//znaleznienie kroku czasowego dt wedlug wzoru (18)

void setU(double u[nx+1][ny+1]){
    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            double x = delta * i;
            double y = delta * j;
            u[i][j] = (1/(2*M_PI*pow(sigma, 2))) * exp(-(pow((x-xA), 2) + pow((y-yA),2))/(2*sigma*sigma));
        }
    }
}//ustalenie warunku poczatkowego wedlug wzoru (17)

void przepiszU(double u1[nx+1][ny+1], double u0[nx+1][ny+1]){
    for(int i=0; i<nx+1; i++){
        for(int j=0; j<ny+1; j++){
            u1[i][j] = u0[i][j];
        }
    }
}

void adwekcja_dyfuzja(double D, double Vx[nx+1][ny+1], double Vy[nx+1][ny+1], double psi[nx+1][ny+1]){

    FILE *fc, *fxsr;
    char filename_c[256];
    char filename_xsr[256];
    sprintf(filename_c, "files/c_%g.txt", D);
    sprintf(filename_xsr, "files/xsr_%g.txt", D);
    fc = fopen(filename_c, "w");
    fxsr = fopen(filename_xsr, "w");
    fclose(fc);
    fclose(fxsr);

    double u0[nx+1][ny+1];
    double u1[nx+1][ny+1];//utworzenie tablic u0 oraz u1

    double dt = find_dt(Vx, Vy);
    setU(u0);

    FILE *gif;
    int t =0;

    printf("rozpoczecie algorytmu AD, D = %g\n", D);
    for(int it = 0; it<IT_MAX; it++){ //glowna petla algorytmu
        przepiszU(u1, u0);

        if(it%1000 == 0) printf("liczenie... %d iteracja\n", it);

        for(int k = 1; k<=20; k++){ // start iteracji Picarda
            for(int i=0; i<=nx; i++){
                for(int j=1; j<=ny-1; j++){
                    if(i>=i1 && i<=i2 && j<=j1){
                        continue;
                    }else if(i==0){
                        u1[i][j] = 1 / (1 + ((2 * D * dt) / (delta * delta))) * (u0[i][j] 
                            - (dt / 2 * Vx[i][j] * ((u0[i + 1][j] - u0[nx][j]) / (2 * delta) + (u1[i + 1][j] - u1[nx][j]) / (2 * delta)))
                            - (dt / 2 * Vy[i][j] * ((u0[i][j + 1] - u0[i][j - 1]) / (2 * delta) + (u1[i][j + 1] - u1[i][j - 1]) / (2 * delta)))
                            + (dt / 2 * D * ((u0[i + 1][j] + u0[nx][j] + u0[i][j + 1] + u0[i][j - 1] - 4 * u0[i][j]) / (delta * delta)
                            + ((u1[i + 1][j] + u1[nx][j] + u1[i][j + 1] + u1[i][j - 1])/ (delta * delta))))); //wzor (9) + periodyczne WB dla i==0 (i-1 zamieniam na nx)
                    }else if(i == nx){
                        u1[i][j] = 1/(1 + ((2 * D *dt) / (delta * delta))) * (u0[i][j] 
                            - (dt / 2 * Vx[i][j] * ((u0[0][j] - u0[i - 1][j]) / (2 * delta) + (u1[0][j] - u1[i - 1][j]) / (2 * delta)))
                            - (dt / 2 * Vy[i][j] * ((u0[i][j + 1] - u0[i][j - 1]) / (2 * delta) + (u1[i][j + 1] - u1[i][j - 1]) / (2 * delta)))
                            + (dt / 2 * D * ((u0[0][j] + u0[i - 1][j] + u0[i][j + 1] + u0[i][j - 1] - 4 * u0[i][j]) / (delta * delta)
                            + ((u1[0][j] + u1[i - 1][j] + u1[i][j + 1] + u1[i][j - 1])/ (delta * delta))))); //wzor (9) + periodyczne WB dla i==nx (i+1 zamieniam na 0)
                    }else{
                        u1[i][j] = 1/(1 + ((2 * D *dt) / (delta * delta))) * (u0[i][j] 
                            - (dt / 2 * Vx[i][j] * ((u0[i + 1][j] - u0[i - 1][j]) / (2 * delta) + (u1[i + 1][j] - u1[i - 1][j]) / (2 * delta))) 
                            - (dt / 2 * Vy[i][j] * ((u0[i][j + 1] - u0[i][j - 1]) / (2 * delta) + (u1[i][j + 1] - u1[i][j - 1]) / (2 * delta)))
                            + (dt / 2 * D * ((u0[i + 1][j] + u0[i - 1][j] + u0[i][j + 1] + u0[i][j - 1] - 4 * u0[i][j]) / (delta * delta)
                            + ((u1[i + 1][j] + u1[i - 1][j] + u1[i][j + 1] + u1[i][j - 1])/ (delta * delta))))); //wzor (9)
                    }

                }
            }
        }

        przepiszU(u0, u1);

        double c = 0, xsr = 0;

        for(int i=0; i<nx+1; i++){
            for(int j=0; j<ny+1; j++){
                c += (u0[i][j] * delta * delta); //obliczenie calki gestosci c
                xsr += ((i * delta) * u0[i][j] * delta * delta); //obliczenie sredniego polozenia pakietu w kierunku x
            }
        }

        fc = fopen(filename_c, "a");
        fxsr = fopen(filename_xsr, "a");
        fprintf(fc, "%lf %lf\n", dt*it, c);
        fprintf(fxsr, "%lf %lf\n", dt*it, xsr);
        fclose(fc);
        fclose(fxsr);

        if(it%(IT_MAX/50) == 0){
            char filename_gif[256];
            sprintf(filename_gif, "files/%gzad5_it=%d.txt", D, t);
            gif = fopen(filename_gif, "w");

            for(int i=0; i<nx+1; i++){
                for(int j=0; j<ny+1; j++){
                    fprintf(gif, "%f\t%f\t%f\n", delta*i, delta*j, u1[i][j]);
                }
                fprintf(gif, "\n");
            }

            fclose(gif);
            t++;
        }//zapis map rozkladu
    }
    printf("zakonczono liczenie\n");
}

int main(){
    mkdir("files", 0777);

    double psi[nx+1][ny+1];

    wczytaj_strumien(psi);

    double Vx[nx+1][ny+1];
    double Vy[nx+1][ny+1];

    setVxy(Vx, Vy, psi);

    adwekcja_dyfuzja(0.0, Vx, Vy, psi); //D = 0 - brak dyfuzji

    adwekcja_dyfuzja(0.1, Vx, Vy, psi); //D = 0.1
}