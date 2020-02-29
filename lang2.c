#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>

// define useful constants
#define PI 3.14159265358979323846
#define KB 1.380658e-23

typedef struct {
    double dt, T, tau, m, gamma, kx, ky, kz;
} params;

typedef struct{
    char valid;
    int sdatadims;
    int slen;
    double * arr;
} data;

double double_rand(){
  double num = rand() / (double) RAND_MAX;
  return num;
}

// Generate a ransom double ~N(0,1) using Box–Muller transform
double random_normal(){
    return sqrt(-2 * log(double_rand())) * cos(2 * PI * double_rand());
}

data make_data(int ddims, int len){
    data s;
    s.valid = 1;
    s.sdatadims = ddims;
    s.slen = len;
    s.arr = (double *)calloc(sizeof(double), len*ddims);
    assert(s.arr != NULL);
    return s;
}

void free_data(data s){
    s.valid = 0;
    free(s.arr);
}

void dump_data(data s, FILE * fp){

    for (int i = 0; i < s.slen; i++){
        for (int j = 0; j < s.sdatadims; j++){
            fprintf(fp, "%.8g\t", *(s.arr + j * s.slen + i));
        }
        fprintf(fp, "\n");
    }
}

void freediff_evolver(double* xi, double * wi, double * time, params p, int num){
    *xi = 0.0f;
    *wi = 0.0f;
    *time = 0.0f;
    double nr;
    double sqdt = sqrt((double)p.dt);
    for (int i = 1; i < num; i++){
        nr = random_normal();
        *(xi + i) = *(xi+i-1) + sqdt * nr;
        *(wi + i) = nr / sqdt;
        *(time + i) = i * p.dt;
    }
}

void brownian_evolver(double * xii, double * xia, double * timetau, params p, int num){
    double D = KB * p.T / p.gamma;
    double gm = 1/p.tau;
    double c1 = (2 + p.dt * gm)/(1 + p.dt * gm);
    double c2 = 1/(1 + p.dt * gm);
    double c3 = (sqrt(2 * KB * p.T * p.gamma) / (p.m * (1 + p.dt * gm))) * pow(p.dt, (double)1.5);
    double c4 = sqrt(2 * D * p.dt);
    // printf("T=%g gamma=%g D=%g c1=%g c2=%g c3=%g c4=%g\n", p.T, p.gamma, D, c1, c2, c3 ,c4);

    *xii = 0.0;
    *xia = 0.0;
    *timetau = 0.0;

    double nr;
    for (int i = 1; i < num; i++){
        nr = random_normal();
        if(i==1) {
            *(xii + i) = *(xia + i - 1) * c1 + c3 * nr ;
        } else {
            *(xii + i) = *(xia + i - 1) * c1 - *(xia + i - 2) * c2 + c3 * nr ;
        }
        *(xia + i) = *(xia + i - 1) + c4 * nr;
        *(timetau + i) = i * p.dt / p.tau;
    }
}

void trap_evolver(double * xi, double * yi, double * zi, double * time, params p, int num){
    *xi = 0.0;
    *yi = 0.0;
    *zi = 0.0;
    *time = 0.0;

    double D = KB * p.T / p.gamma;
    double c = sqrt(2 * D * p.dt);

    for (int i = 1; i < num; i++)
    {
        *(xi + i) = *(xi + i - 1) * (1 - p.kx * p.dt/p.gamma) + c * random_normal();
        *(yi + i) = *(yi + i - 1) * (1 - p.ky * p.dt/p.gamma) + c * random_normal();
        *(zi + i) = *(zi + i - 1) * (1 - p.kz * p.dt/p.gamma) + c * random_normal();
        *(time + i) = p.dt * i;
    }
}
    
int main(int argc, char *argv[]) {
    srand(time(0));

    // parse params
    params p;

    p.dt = 0.1f;
    if (argc >= 3) sscanf(argv[2], "%lg", &p.dt);
    
    int nsteps = 500;
    if (argc >= 4) sscanf(argv[3], "%d", &nsteps);
    

    char fn[100];

    if (0 == strcmp(argv[1], "freediff")){
        printf("Running free diffusion with params nsteps=%d - dt=%g\n", nsteps, p.dt);
        sprintf(fn, "FreeDiff_%i_%.3f.out", nsteps, p.dt);
        FILE *fp = fopen(fn,"w");
        data d = make_data(3, nsteps);
        freediff_evolver(d.arr, d.arr + nsteps, d.arr + 2 * nsteps, p, nsteps);
        dump_data(d, fp);
        fclose(fp);
        free_data(d);
    } else if(0 == strcmp(argv[1], "browndiff")){
        //// System params
        double eta = 0.001;
        double R=1E-6;
        p.gamma = 6 * PI * eta * R;
        p.m = 1.1e-14;
        p.tau = p.m/p.gamma;
        p.T = 300;

        printf("Running brownian diffusion with params nsteps=%d - dt=%g [tau=%g]\n", nsteps, p.dt, p.tau);
        sprintf(fn, "BrownDiff_%i_%.3g.out", nsteps, p.dt);
        FILE *fp = fopen(fn,"w");
        data d = make_data(3, nsteps);
        brownian_evolver(d.arr, d.arr + nsteps, d.arr + 2 * nsteps, p, nsteps);
        dump_data(d, fp);
        fclose(fp);
        free_data(d);
    }else if(0 == strcmp(argv[1], "opttrap")){
        //// System params
        double eta = 0.001;
        double R=1E-6;
        p.T = 300;
        p.gamma = 6 * PI * eta * R;
        p.m = 1.1E-14;
        p.tau = p.m/p.gamma;

        p.kx = 1E-6;
        p.ky = 1E-6;
        p.kz = 1E-6;
        if (argc >= 5) sscanf(argv[4], "%lg_%lg_%lg", &p.kx, &p.ky, &p.kz);

        printf("Running optical trap with params n=%d - dt=%g - [kx=%g, ky=%g, kz=%g]\n", nsteps, p.dt, p.kx, p.ky, p.kz);
        sprintf(fn, "OptTrap_%i_%.3g_%.3g_[%.3g_%.3g_%.3g].out", nsteps, p.dt, p.tau, p.kx, p.ky, p.kz);
        FILE *fp = fopen(fn,"w");
        data d = make_data(4, nsteps);
        trap_evolver(d.arr, d.arr + nsteps, d.arr + 2 * nsteps, d.arr + 3 * nsteps, p, nsteps);
        dump_data(d, fp);
        fclose(fp);
        free_data(d);
    } else {
        printf("use : lang [freediff / browndiff / opttrap] [dT(=0.1)] [nsteps (=500)] [kx_ky_kz]\n");
    }
    return 0;
}