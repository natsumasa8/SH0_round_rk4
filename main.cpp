#include <iostream>
#include <Eigen/Core>
#include <cmath>
#include "Calcdif.h"
#include "ioput.h"
#include "SHsimulate.h"
#include "initial.h"

int main(){
    Calcdif clc;
    clc.readinput();

    SHsimulate shs;
    shs.ConstractParameters();

    const int Nx = clc.Nx;
    IOput ioput;
    ioput.Nx = Nx;
    std::cout << "Nx in main.cpp is " << Nx << std::endl;

    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(Nx,Nx);
    Eigen::MatrixXd u_theta = Eigen::MatrixXd::Zero(Nx,Nx);
    Eigen::MatrixXd u_grad_x = Eigen::MatrixXd::Zero(Nx,Nx);
    Eigen::MatrixXd u_grad_y = Eigen::MatrixXd::Zero(Nx,Nx);
    Eigen::MatrixXd u_old = Eigen::MatrixXd::Zero(Nx,Nx);

    Initial ini;
    ini.set_type1(u_old, clc.initial_amp); // type1: random

    shs.round_boundary(u_old, clc.u_outside); //setting round boundary condition

    for(int n=1;n<=clc.Nt;n++){

        u = shs.CalcSH_rk4(u_old,clc);

        shs.round_boundary(u, clc.u_outside); //setting round boundary condition
        
        if(n%clc.N_output==0 || n==1){
            //calculation of gradient
            u_grad_x = clc.grad_x(u);
            u_grad_y = clc.grad_y(u);
            u_theta = clc.theta(u_grad_x, u_grad_y);
            //calculation of energy(density)
            Eigen::MatrixXd energy_density = shs.CalcED(u, u_grad_x, u_grad_y, clc);                
            double energy = shs.CalcE(energy_density, clc);
            Eigen::VectorXd energy_term = shs.CalcED_term(u, u_grad_x, u_grad_y, clc);

            //output ->vtk, txt
            char filename[100];
            sprintf(filename, "d%08d.vtk",n);
            ioput.output(u, u_theta, energy_density, filename);
            shs.output_energy(n, energy, energy_term, "energy.txt", clc);

            std::cout << "step: " << n << ", enrgy: " << energy << std::endl;
        }
        u_old = u;
    }
    std::cout << "done" << std::endl;
    return 0;
}