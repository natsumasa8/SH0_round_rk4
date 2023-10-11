#include <Eigen/Core>
#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include "Calcdif.h"

// functions
void Calcdif::readinput(){
	std::cout << "Put input file below." << std::endl;
	std::cin >> FilePath;
	std::string fname = FilePath;
	std::ifstream fin(fname);
	if(!fin){
		std::cout << "Error: cannot open input file." << std::endl;
		exit(1);
	}
//reading inputfile
	std::string dummy, check;
	fin >> dummy; fin >> dt; std::cout << "dt: " << dt << std::endl;
	fin >> dummy; fin >> Nt; std::cout << "Nt: " << Nt << std::endl;
	fin >> dummy; fin >> N_output; std::cout << "N_output: " << N_output << std::endl;
	
	fin >> dummy; fin >> dx; std::cout << "dx: " << dx << std::endl;
	fin >> dummy; fin >> Nx; std::cout << "Nx: " << Nx << std::endl;
	
	fin >> dummy; fin >> u_outside; std::cout << "u(outside): " << u_outside << std::endl;

	fin >> check;
	if(check != "End"){
		std::cout << "Error: invalid input file." << std::endl;
		exit(1);
	}
	fin.close();
	std::cout << "Input has successfuly done." << std::endl;

//stability check
	double stab = 2 * dt / std::pow(dx, 4);
	if(stab > 0.5){
		std::cout << "Warning : Stability condition is not satisfied." << std::endl;
		std::cout << "continue?[y/n]" << std::endl;
		fin >> dummy;
		if(dummy=="n"){
			exit(1);
		}
	}
}

double Calcdif::Dx2(int p, int q, Eigen::MatrixXd& s){
    double tmp = s((p+1)%Nx, q) + s((p-1+Nx)%Nx, q) - 2.0 * s(p, q); 
    tmp /= dx * dx;
    return tmp;
}

double Calcdif::Dy2(int p, int q, Eigen::MatrixXd& s){
    double tmp = s(p, (q+1)%Nx) + s(p, (q-1+Nx)%Nx) - 2.0 * s(p, q);
    tmp /= dx * dx;
    return tmp;
}

double Calcdif::Dx4(int p, int q, Eigen::MatrixXd& s){
    double tmp = s((p+2)%Nx, q) - 4.0 * s((p+1)%Nx, q) + 6.0 * s(p, q) - 4.0 *  s((p-1+Nx)%Nx, q) + s((p-2+Nx)%Nx, q);
    tmp /= dx * dx * dx * dx;
    return tmp;
}

double Calcdif::Dy4(int p, int q, Eigen::MatrixXd& s){
    double tmp = s(p, (q+2)%Nx) - 4.0 * s(p, (q+1)%Nx) + 6.0 * s(p, q) - 4.0 *  s(p, (q-1+Nx)%Nx) + s(p, (q-2+Nx)%Nx);
    tmp /= dx * dx * dx * dx;
    return tmp;
}

double Calcdif::Dx2y2(int p, int q, Eigen::MatrixXd& s){
    double tmp = 0.0;
    tmp += s((p+1)%Nx, (q+1)%Nx);
	tmp -= 2.0 * s((p+1)%Nx, q);
    tmp += s((p+1)%Nx, (q-1+Nx)%Nx);
	tmp -= 2.0 * s(p, (q+1)%Nx);
	tmp += 4.0 * s(p, q);
    tmp -= 2.0 * s(p, (q-1+Nx)%Nx);
	tmp += s((p-1+Nx)%Nx, (q+1)%Nx);
	tmp -= 2.0 * s((p-1+Nx)%Nx, q);
    tmp += s((p-1+Nx)%Nx, (q-1+Nx)%Nx);
    
    tmp /= dx * dx * dx * dx;
	return tmp;
}

double Calcdif::laplacian(int p, int q, Eigen::MatrixXd& s){
    return Dx2(p,q,s) + Dy2(p,q,s);
}

double Calcdif::laplacian2(int p, int q, Eigen::MatrixXd& s){
    return Dx4(p,q,s) + 2.0 * Dx2y2(p,q,s) + Dy4(p,q,s);
}

//calculated gradient.
Eigen::MatrixXd Calcdif::grad_x(Eigen::MatrixXd& s){
	Eigen::MatrixXd g_x = Eigen::MatrixXd::Zero(Nx, Nx);
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Nx;j++){
			double tmp = s((i+1)%Nx,j) - s((i-1+Nx)%Nx,j);
			tmp /= 2.0 * dx;
			g_x(i,j) = tmp;
		}
	}
	return g_x;
}

Eigen::MatrixXd Calcdif::grad_y(Eigen::MatrixXd& s){
	Eigen::MatrixXd g_y = Eigen::MatrixXd::Zero(Nx, Nx);
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Nx;j++){
			double tmp = s(i,(j+1)%Nx) - s(i,(j-1+Nx)%Nx);
			tmp /= 2.0 * dx;
			g_y(i,j) = tmp;
		}
	}
	return g_y;
}

Eigen::MatrixXd Calcdif::theta(Eigen::MatrixXd& gx, Eigen::MatrixXd& gy){
	Eigen::MatrixXd tht = Eigen::MatrixXd::Zero(Nx,Nx);
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Nx;j++){
			double tmp = 0.0;
			tmp = atan2(gy(i,j),gx(i,j));
			if(tmp<0){
				tht(i,j) = tmp + M_PI;
			}else{
				tht(i,j);
			}
		}
	}
	return tht;
}