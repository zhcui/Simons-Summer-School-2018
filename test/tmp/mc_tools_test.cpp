/******
2D Ising model C tools, for Simons Summer School 2018
Author: Zhihao Cui
        Gaurav Harsha 
******/

#include <stdio.h>
#include <iostream>
#include <random>
#include <math.h>

void c_mc_step(int* A, const int & length, const double* pflip, double & energy, double & mag){


    int N = length * length;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis_int(0, length - 1);
    std::uniform_real_distribution<double> dis_real(0.0, 1.0);

    for(int i = 0; i < N; ++i){
        
        // randomly select a site
        int x = dis_int(gen);
        int y = dis_int(gen);
        
        // local interaction count, PBC is used
        int env_factor = A[x * length + y] * (A[((x - 1 + length)%length) * length + y] + A[((x + 1)%length) * length + y] + 
                            A[x * length + (y - 1 + length)%length] + A[x * length + (y + 1)%length]);
        // flip
        double r = dis_real(gen);
        if(r <= pflip[env_factor + 4]){
            A[x * length + y] *= -1;
            energy += env_factor * 2.0;
            mag += A[x * length + y] * 2.0;
        }
    }
}

int main(){
    using namespace std;

    int length = 10;
    double temp = 6.0, energy = 0.0, mag = 0.0;
    double pflip[9] = {0.0};
    int A[length * length] = {0}; 
    for (int i = 0; i < length; ++i){
        for (int j = 0; j < length; ++j){
            A[i * length + j] = 1;
        }
    }
    int index[5] = {-4, -2, 0, 2, 4};
    for (int i = 0; i < 5; ++i){
        int idx = index[i];
        pflip[idx + 4] = exp(-idx * 2.0 / float(temp));
    }
    c_mc_step(A, length, pflip, energy, mag);
    for (int i = 0; i < length; ++i){
        for (int j = 0; j < length; ++j){
            printf("%2d ", A[i * length + j]);
        }
        cout<<"\n";
    }
    printf("%f", energy);
    printf("%f", mag);
        

    return 0;
}
