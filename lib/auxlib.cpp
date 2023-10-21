#include <stdlib.h>
#include <stdio.h>

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