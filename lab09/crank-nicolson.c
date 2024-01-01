#include<stdio.h>
#include<math.h>

#include<sys/stat.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define nx 40
#define ny 40
#define N ((nx+1)*(ny+1))
#define delta 1.0
#define dt 1.0
#define T_A 40
#define T_B 0
#define T_C 30
#define T_D 0
#define k_B 0.1
#define k_D 0.6
#define IT_MAX 2000

void fill_matrices(gsl_matrix *A, gsl_matrix *B, gsl_vector *C){
    for(int i=1; i<nx; i++){
        for(int j=1; j<ny; j++){
            int l = i + j * (nx+1);
            double val = dt/(2.0*delta*delta);

            gsl_matrix_set(A, l, l-nx-1, val);
            gsl_matrix_set(A, l, l-1, val);
            gsl_matrix_set(A, l, l+1, val);
            gsl_matrix_set(A, l, l+nx+1, val);

            val = -(2.0*dt)/(delta*delta) - 1;

            gsl_matrix_set(A, l, l, val);

            val = -dt/(2.0*delta*delta);

            gsl_matrix_set(B, l, l-nx-1, val);
            gsl_matrix_set(B, l, l-1, val);
            gsl_matrix_set(B, l, l+1, val);
            gsl_matrix_set(B, l, l+nx+1, val);

            val = (2.0*dt)/(delta*delta) - 1;

            gsl_matrix_set(B, l, l, val);
        }
    }

    for(int i=0; i<=nx; i+=nx){
        for(int j=0; j<=ny; j++){
            int l = i + j * (nx+1);
            gsl_matrix_set(A, l, l, 1.0);
            gsl_matrix_set(B, l, l, 1.0);
            gsl_vector_set(C, l, 0.0);
        }
    }

    for(int i=1; i<=nx-1; i++){
        int j = ny;
        int l = i + j * (nx+1);
        double val = -1.0/(k_B*delta);

        gsl_matrix_set(A, l, l-nx-1, val);

        val = 1.0 + 1.0/(k_B*delta);

        gsl_matrix_set(A, l, l, val);
        gsl_vector_set(C, l, T_B);
        for(int k=0; k<N; k++){
            gsl_matrix_set(B, l, k, 0.0);
        }
    }

    for(int i=1; i<=nx-1; i++){
        int j = 0;
        int l = i + j * (nx+1);
        double val = 1.0 + 1.0/(k_D*delta);

        gsl_matrix_set(A, l, l, val);

        val = -1.0/(k_D*delta);

        gsl_matrix_set(A, l, l+nx+1, val);
        gsl_vector_set(C, l, T_D);
        for(int k=0; k<N; k++){
            gsl_matrix_set(B, l, k, 0.0);
        }
    }
}

void fill_vector_T(gsl_vector *T){
    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++){
            int l = i + j * (nx+1);

            gsl_vector_set(T, l, 0);

            if(i==0){
                gsl_vector_set(T, l, T_A);
            }

            if(i==nx){
                gsl_vector_set(T, l, T_C);
            }
        }
    }
}

void save_map(gsl_vector *T, int it){
    FILE *f_T, *f_NT;
    char filename_T[50];
    char filename_NT[50];
    sprintf(filename_T, "files/T_%d.txt", it);
    sprintf(filename_NT, "files/NT_%d.txt", it);

    f_T = fopen(filename_T, "w");
    f_NT = fopen(filename_NT, "w");

    for(int i=1; i<nx; i++){
        for(int j=1; j<ny; j++){
            int l = i + j * (nx+1);
            int x = delta * i;
            int y = delta * j;
            fprintf(f_T, "%d %d %lf\n", x, y, gsl_vector_get(T, l));

            double nabla_T;
            nabla_T = ((gsl_vector_get(T, l + 1) - 2.0*gsl_vector_get(T, l) + gsl_vector_get(T, l-1))/pow(delta, 2)) +
            ((gsl_vector_get(T, l + nx + 1) - 2.0*gsl_vector_get(T, l) + gsl_vector_get(T, l- nx -1) )/pow(delta, 2));

            fprintf(f_NT, "%d %d %.10lf\n", x, y, nabla_T);
        }
        fprintf(f_T, "\n");
        fprintf(f_NT, "\n");
    }

    fclose(f_NT);
    fclose(f_T);
}

void CN_algorithm(gsl_matrix *A, gsl_matrix *B, gsl_vector *c, gsl_vector *T_0){
    int signum = 0;
    gsl_permutation *p = gsl_permutation_calloc(N);

    gsl_vector *T_n = gsl_vector_calloc(N);
    gsl_vector *T_n1 = gsl_vector_calloc(N);
    gsl_vector_memcpy(T_n, T_0);
    gsl_vector_memcpy(T_n1, T_0);

    gsl_linalg_LU_decomp(A, p, &signum);

    for(int it=0; it<=IT_MAX; it++){
        gsl_vector *d = gsl_vector_calloc(N);
        gsl_blas_dgemv(CblasNoTrans, 1.0, B, T_n, 0.0, d);
        gsl_vector_add(d, c);

        gsl_linalg_LU_solve(A, p, d, T_n1);

        gsl_vector_memcpy(T_n, T_n1);

        gsl_vector_free(d);

        if(it==100 || it==200 || it==500 || it==1000 || it==2000){
            save_map(T_n, it);
        }

        if(it%100==0) printf("liczenie... %d iteracja\n", it);
    }

    gsl_permutation_free(p);
    gsl_vector_free(T_n);
    gsl_vector_free(T_n1);
}

int main(){
    mkdir("files", 0777);

    gsl_matrix *A = gsl_matrix_calloc(N, N);
    gsl_matrix *B = gsl_matrix_calloc(N, N);
    gsl_vector *c = gsl_vector_calloc(N);

    fill_matrices(A, B, c);

    gsl_vector *T_0 = gsl_vector_calloc(N);

    fill_vector_T(T_0);

    CN_algorithm(A, B, c, T_0);

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(c);
    gsl_vector_free(T_0);

    mkdir("img", 0777);

    return 0;
}