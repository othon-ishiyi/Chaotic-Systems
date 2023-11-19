#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.14159

void henon_map(float a, float b, float &x, float &y) { 
    float x_aux = a - x*x + b*y;
    y = x;
    x = x_aux;
}

extern "C"{
    float* henon_attractor(float a, float b, float x0, float y0, float xmin, float xmax, float ymin, float ymax, int niter) {
        float x = x0;
        float y = y0;

        float *attractor = (float *) malloc(2*niter*sizeof(float));
        int i = 0;

        while(i < niter){
            henon_map(a, b, x, y);
            if ((x > xmin) && (x < xmax) && (y > ymin) && (y < ymax)){
                attractor[i] = x;
                attractor[i+niter] = y;
                i++;
            }
        }

        return attractor;
    }
}

extern "C"{
    void free_vector(float* v) {
        free(v);
    } 
}

double dtheta2_dt(double t, double theta, double dtheta) {
    return -0.2*dtheta - sin(theta) + 2*cos(t);
}

double* solve_pendulum(double theta0, double dtheta0, double dt, double t_max) {
    int size = ceil(t_max/dt)+1;
    double t_ant, t;
    double k1_theta, k2_theta, k3_theta, k4_theta;
    double k1_dtheta, k2_dtheta, k3_dtheta, k4_dtheta;
    double theta, dtheta;

    theta = theta0;
    dtheta = dtheta0;
    t_ant = 0;

    for(int i=0; i<size-1; i++) {
        t = t_ant + dt;
        k1_theta = dt*dtheta;
        k1_dtheta = dt*dtheta2_dt(t_ant, theta, dtheta);
        k2_theta = dt*(dtheta+k1_dtheta/2);
        k2_dtheta = dt*dtheta2_dt(t_ant+dt/2, theta+k1_theta/2, dtheta+k1_dtheta/2);
        k3_theta = dt*(dtheta+k2_dtheta/2);
        k3_dtheta = dt*dtheta2_dt(t_ant+dt/2, theta+k2_theta/2, dtheta+k2_dtheta/2);
        k4_theta = dt*(dtheta+k3_dtheta);
        k4_dtheta = dt*dtheta2_dt(t, theta+k3_theta, dtheta+k3_dtheta);
        theta += (k1_theta+2*k2_theta+2*k3_theta+k4_theta)/6;
        dtheta += (k1_dtheta+2*k2_dtheta+2*k3_dtheta+k4_dtheta)/6;
        t_ant = t;
    }

    theta = fmod(theta, 2*PI);
    if(theta > PI) {
        theta = theta - 2*PI;
    }
    double* sol = (double *) malloc(2*sizeof(double));
    sol[0] = theta;
    sol[1] = dtheta;

    return sol;
}

double* solve_pendulum_euler(double theta0, double dtheta0, double dt, double t_max) {
    int size = ceil(t_max/dt)+1;
    double t_ant, t;
    
    double *theta = (double *) malloc(size*sizeof(double));
    double *dtheta = (double *) malloc(size*sizeof(double));

    theta[0] = theta0;
    dtheta[0] = dtheta0;
    t_ant = 0;

    for(int i = 0; i < size-1; i++) {
        t = t_ant + dt;
        theta[i+1] = theta[i] + dt*dtheta[i];
        dtheta[i+1] = dtheta[i] + dt*dtheta2_dt(t_ant, theta[i], dtheta[i]);
        t_ant = t;
    }

    double* sol = (double *) malloc(2*sizeof(double));
    sol[0] = theta[size-1];
    sol[0] = fmod(sol[0], 2*PI);
    if(sol[0] > PI) {
        sol[0] = sol[0] - 2*PI;
    }
    sol[1] = dtheta[size-1];
    free(theta);
    free(dtheta);

    return sol;

}

extern "C"{
    int classif_by_attractor(double theta0, double dtheta0) {
        double* sol = solve_pendulum(theta0, dtheta0, 0.01, 100*PI);
        if (((sol[0]+0.477)*(sol[0]+0.477) + (sol[1]+0.609)*(sol[1]+0.609)) < 0.01) {
            return 0;
        }
        else if (((sol[0]+0.471)*(sol[0]+0.471) + (sol[1]-2.037)*(sol[1]-2.037)) < 0.01) {
            return 1;
        }
        else{
            return -1;
        }
    }
}

extern "C"{
    int classif_by_attractor_light(double theta0, double dtheta0) {
        double* sol = solve_pendulum(theta0, dtheta0, 0.01, 30*PI);
        if (sol[1] < 0.5) {
            return 0;
        }
        else{
            return 1;
        }
    }
}



int main() {
    int a, b;
    for(float i = 0; i <= 100; i++) {
        for(float j = 0; j <= 100; j++) {
            a = classif_by_attractor(-3.0 + i*3.0/50.0, -3.0 + j*3.0/50.0);
            b = classif_by_attractor_light(-3.0 + i*3.0/50.0, -3.0 + j*3.0/50.0);
            if(a != b){
                printf("%d ", a);
            }
        }
        printf("\n");
    }
    return 0;
}


