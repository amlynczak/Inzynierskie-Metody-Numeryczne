#include<stdio.h>
#include<math.h>
#include<sys/stat.h>

/*ustalenie wartosci parametrow*/
#define nx 150
#define nt 1000
#define delta 0.1
#define dt 0.05
#define xA 7.5
#define sigma 0.5

#define t_max (dt*nt)
#define xF 2.5
/**/

double delta_kroneker(double x){
    return fabs(x - xF) < 1e-8 ? 1.0 : 0.0; //delta Kroneckera
}

double a_F(double x, double t){
    return cos((50*t)/t_max)*delta_kroneker(x); //wzor (3)
}

void warunki_brzeg(double *u, double *v){ //sztywne warunki brzegowe, wzory (12) oraz (13)
    u[0] = 0;
    u[nx] = 0;
    v[0] = 0;
    v[nx] = 0;
}

void warunki_start(double *u, double *v){ //warunki poczatkowe, wzory (14) i (15)
    double x;
    for(int i=1; i<nx; i++){
        x = i*delta;
        u[i] = exp (- (pow(x - xA, 2) / (2.0* pow(sigma,2))));
        v[i] = 0;
    }
}

void set_a(double *a, double *u, double *u0, double alpha, double beta, double t){ //ustalenie wartosci w tablicy a, wzor (11)
    for(int i=1; i<nx; i++){
        double x = i*delta;
        a[i] = (u[i+1] - 2.0 * u[i] + u[i-1]) / (pow(delta,2)) - beta * ((u[i] - u0[i]) / dt) + alpha * a_F(x, t);
    }
}

double energia(double *u, double *v){ //obliczenie wartosci energii, wzor (17)
    double var1, sum = 0;

    var1 = (delta / 4.0) * (pow((u[1] - u[0])/delta, 2) + pow((u[nx] - u[nx-1])/delta, 2));
    
    for(int i=1; i<=nx-1; i++){
        sum = sum + delta/2.0 * (pow(v[i],2) + (pow((u[i+1] - u[i-1])/(2.0*delta), 2)));
    }

    return var1 + sum;
}

void rownianie_falowe(double alpha, double beta){ //algorytm rowania falowego
    char filename1[50];
    sprintf(filename1, "files/E_a%.01f_b%.01f.txt", alpha, beta);
    FILE *file1;
    file1 = fopen(filename1, "w");

    char filename2[50];
    sprintf(filename2, "files/map_a%.01f_b%.01f.txt", alpha, beta);
    FILE *file2;
    file2 = fopen(filename2, "w");

    double u0[nx+1];
    double u[nx+1];
    double v[nx+1];
    double vp[nx+1];
    double a[nx+1]; //inicjalizacja tablic 1D
    double E = 0;

    warunki_brzeg(u, v);

    if(beta != 1.0 || alpha != 1.0) { //dla ostatniego punktu przyjmujemy inne warunki poczatkowe
        warunki_start(u, v);
    }else{
        for(int i=1; i<nx; i++){
            u[i] = 0;
            v[i] = 0;
        }
    }

    for(int i=0; i<=nx; i++){
        u0[i] = u[i];
    } //zachowanie poprzedniego wyniku

    set_a(a, u, u0, alpha, beta, 0); //inicjalizacja a

    for(int n = 1; n <= nt; n++){ //glowna petla algorytmu
        double t = dt*n;
        for(int i=1; i<nx; i++){
            vp[i] = v[i] + (dt/2.0)*a[i];
            u0[i] = u[i];
            u[i] = u[i] + dt * vp[i];
        }

        for(int i=1; i<nx; i++){
            set_a(a, u, u0, alpha, beta, t);
            v[i] = vp[i] + (dt/2.0)*a[i];
        }

        E = energia(u, v);

        fprintf(file1, "%g %g \n", t, E);
        for(int i=0; i<=nx; i++){
            fprintf(file2, "%lf %lf %g\n", t, i*delta, u[i]); //zapis do pliku
        }
    }

    fclose(file1);
    fclose(file2);
}

int main(){
    mkdir("files", 0777);

    double alpha, beta;

    ////
    alpha = 0.0;
    beta = 0.0;
    rownianie_falowe(alpha, beta);

    ////
    alpha = 0.0;
    beta = 0.1;
    rownianie_falowe(alpha, beta);

    ////
    alpha = 0.0;
    beta = 1.0;
    rownianie_falowe(alpha, beta);

    ////
    alpha = 1.0;
    beta = 1.0;
    rownianie_falowe(alpha, beta);

}