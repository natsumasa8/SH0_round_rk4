#include <Eigen/Core>
#include <cmath>
#include <iostream>
#include <fstream>
#include "SHsimulate.h"
#include "Calcdif.h"

void SHsimulate::ConstractParameters(){
    std::cout << "Construct SH-parameters...put filepath below." << std::endl;
    std::cin >> FilePath;

    std::string fname = FilePath;
    std::ifstream fin(fname);
    if(!fin){
        std::cout << "Error: cannot open parameter file." << std::endl;
        exit(1);
    }

    std::string dummy;
    fin >> dummy; fin >> a; std::cout << "a: " << a << std::endl;
    fin >> dummy; fin >> b; std::cout << "b: " << b << std::endl;
    fin >> dummy; fin >> c; std::cout << "c: " << c << std::endl;
    
    std::string check;

    fin >> check;
    if(check != "End"){
        std::cout << "Error: invalid parameter file." << std::endl;
        exit(1);
    }
    fin.close();
}

//note: clc -> s.row() etc...
Eigen::MatrixXd SHsimulate::CalcSH(Eigen::MatrixXd& s, Calcdif& clc){
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(clc.Nx,clc.Nx);
#pragma omp parallel for
    for (int p=0;p<clc.Nx;p++){
        for (int q=0;q<clc.Nx;q++){
            tmp(p,q) += -clc.laplacian(p,q,s);
            tmp(p,q) += -clc.laplacian2(p,q,s);
            tmp(p,q) += -a * s(p,q);
            tmp(p,q) += -b * std::pow(s(p,q),2.0);
            tmp(p,q) += -c * std::pow(s(p,q),3.0);
            
            tmp(p,q) *= clc.dt;
            tmp(p,q) += s(p,q);
            // tmp(p,q) += rand; //gaussian white noise ?? every time? 
        }
    }
    return tmp;
}

//初期緩和計算用
//note: clc -> s.row() etc...
Eigen::MatrixXd SHsimulate::CalcRX(Eigen::MatrixXd& s, Calcdif& clc){
    Eigen::MatrixXd rx = Eigen::MatrixXd::Zero(clc.Nx,clc.Nx);
#pragma omp parallel for
    for (int p=0;p<clc.Nx;p++){
        for (int q=0;q<clc.Nx;q++){
            rx(p,q) += -clc.laplacian(p,q,s);
            rx(p,q) += -clc.laplacian2(p,q,s);
            rx(p,q) += -a * s(p,q);
            rx(p,q) += -b * std::pow(s(p,q),2.0);
            rx(p,q) += -c * std::pow(s(p,q),3.0);
            
            rx(p,q) *= clc.dt;
            rx(p,q) += s(p,q);
        }
    }
    return rx;
}

//round boundary condition
void SHsimulate::round_boundary(Eigen::MatrixXd& s, double b){
    // double b is outside value(constant).
    double r = (double)s.rows() / 2.0 - 2.0;
#pragma omp parallel for	
    for(int i=0;i<s.rows();i++){
        for(int j=0;j<s.rows();j++){
            // メッシュ座標から円の中心までの距離を計算
            double distance = std::pow(i - s.rows()/2, 2) + std::pow(j - s.rows()/2, 2);
            if (r*r <= distance) {
                s(i,j) = b;
            }
        }
    }   
}

//エネルギー密度
Eigen::MatrixXd SHsimulate::CalcED(Eigen::MatrixXd& s, Eigen::MatrixXd& gs_x, Eigen::MatrixXd& gs_y, Calcdif& clc){
    Eigen::MatrixXd tmp = Eigen::MatrixXd::Zero(clc.Nx,clc.Nx);
#pragma omp paralle for
	for (int i=0;i<clc.Nx;i++){
        for (int j=0;j<clc.Nx;j++){
            tmp(i,j) += c * std::pow(s(i,j),4.0) / 4.0;
            tmp(i,j) += b * std::pow(s(i,j),3.0) / 3.0;
            tmp(i,j) += a * std::pow(s(i,j),2.0) / 2.0;
			tmp(i,j) -= (gs_x(i,j)*gs_x(i,j) + gs_y(i,j)*gs_y(i,j)) / 2.0; //この/2がダメなのか？？
			tmp(i,j) += std::pow(clc.laplacian(i,j,s),2.0) / 2.0;
        }        
    }
	return tmp;
}

//エネルギー密度（各項）
Eigen::VectorXd SHsimulate::CalcED(Eigen::MatrixXd& s, Eigen::MatrixXd& gs_x, Eigen::MatrixXd& gs_y){
    Eigen::VectorXd ene_every= Eigen::Vector::Zero(5);
    #pragma parallel for 
    for (int i=0;i<clc.Nx;i++){
        for (int j=0;j<clc.Nx;j++){
            
        }
    }
}

//エネルギー積分（Riemann積分）
double SHsimulate::CalcE(Eigen::MatrixXd& s, Calcdif& clc){
	double sum = 0.0;
	for (int i=0;i<clc.Nx;i++){
		for (int j=0;j<clc.Nx;j++){
			sum += s(i,j);
		}
	}
	sum *= std::pow(clc.dx, 2);
	return sum;
}

void SHsimulate::output_energy(int step, double energy, const char* filename, Calcdif& clc){
	std::ofstream ofs(filename, std::ios::app);
	if(!ofs){
		std::cout << "ioput.cpp:output_energy error: unable to open " << filename << "." << std::endl;
	}
	if (step==1){
		ofs << "dt: " << clc.dt << std::endl;
		ofs << "Nt: " << clc.Nt << std::endl;
		ofs << "N_output: " << clc.N_output << std::endl;
		ofs << "dx: " << clc.dx << std::endl;
		ofs << "Nx: " << clc.Nx << std::endl;
        ofs << "initial_amp: " << clc.initial_amp << std::endl; 
		ofs << "STEP,ENERGY" << std::endl;
	}
	ofs << step << "," << energy << std::endl;
	ofs.close();
}