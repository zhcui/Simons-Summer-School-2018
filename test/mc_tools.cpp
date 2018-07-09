/******
2D Ising model C tools, for Simons Summer School 2018
Author: Zhihao Cui
        Gaurav Harsha 
******/

#include <random>
#include <stdio.h>

void c_mc_step(int* A, const int & length, const double* pflip, double & energy, double & mag){


    int N = length * length;
    int x, y, env_factor;
    double r;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis_int(0, length - 1);
    std::uniform_real_distribution<double> dis_real(0.0, 1.0);


    for(int i = 0; i < N; ++i){
        
        // randomly select a site
        x = dis_int(gen);
        y = dis_int(gen);
        
        // local interaction count, PBC is used
        env_factor = A[x * length + y] * (A[((x - 1 + length)%length) * length + y] + A[((x + 1)%length) * length + y] + 
                            A[x * length + (y - 1 + length)%length] + A[x * length + (y + 1)%length]);
        // flip
        r = dis_real(gen);
        if(r <= pflip[(env_factor+9)%9]){
            A[x * length + y] *= -1;
            energy += (env_factor * 2.0);
            mag += (A[x * length + y] * 2.0);
        }
    }
}

void c_mc_step_pre(int* A, const int & length, const double* pflip){


    int N = length * length;
    int x, y, env_factor;
    double r;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis_int(0, length - 1);
    std::uniform_real_distribution<double> dis_real(0.0, 1.0);

    for(int i = 0; i < N; ++i){
        
        // randomly select a site
        x = dis_int(gen);
        y = dis_int(gen);
        
        // local interaction count, PBC is used
        env_factor = A[x * length + y] * (A[((x - 1 + length)%length) * length + y] + A[((x + 1)%length) * length + y] + 
                            A[x * length + (y - 1 + length)%length] + A[x * length + (y + 1)%length]);
        // flip
        r = dis_real(gen);
        if(r <= pflip[(env_factor+9)%9]){
            A[x * length + y] *= -1;
        }
    }
}
