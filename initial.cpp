#include <iostream>
#include <random>
#include <Eigen/Core>
#include <fstream>
#include "initial.h"

void Initial::readinitial(){
	std::cout << "Put initial condition file below." << std::endl;

	std::string FilePath;
	std::cin >> FilePath;
	std::ifstream fin(FilePath);
	if(!fin){
		std::cout << "Error: cannot open input file." << std::endl;
		exit(1);
	}

//reading inputfile
	std::string dummy, check;
	fin >> dummy; fin >> type; std::cout << "type(initial condition): " << type << std::endl;
	if(type=="type1"){
		//no processing
	}else if(type=="type2"){
		for(int i=0;i<defect_num;i++){
		fin >> dummy; fin >> X[i] >> Y[i]; std::cout << "(X1,Y1): " << X[i] << " " << Y[i] << std::endl;
		}
		fin >> dummy; fin >> r; std::cout << "raduis: " << r << std::endl;
	}else{
		std::cout << "Error: invalid input of 'type'." << std::endl;
	}
	fin >> dummy; fin >> amplitude; std::cout << "amplitude: " << amplitude << std::endl;

	fin >> check;
	if(check != "End"){
		std::cout << "Error: invalid input file." << std::endl;
		exit(1);
	}
	fin.close();
	std::cout << "Input has successfuly done." << std::endl;
}

//random-------------------------------------------------------
void Initial::set_type1(Eigen::MatrixXd& s, double a){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> rand001(-1.0,1.0);

#pragma omp parallel for
    for(int i=0;i<s.rows();i++){
		for(int j=0;j<s.rows();j++){
			s(i,j) += a * rand001(mt);
		}
	}
}
//vortex-------------------------------------------------------
void Initial::set_type2(Eigen::MatrixXd& s, int X, int Y){
	//note : r is pixels not real distance
#pragma omp parallel for
	for(int i=0;i<s.rows();i++){
		for(int j=0;j<s.rows();j++){
			double distance = std::pow(i-X,2) + std::pow(j-Y,2);
			if (((r-1)*(r-1))<= distance && distance<=((r+1)*(r+1))) {
                s(i,j) = amplitude;
			}
		}
	}
}

//後で
void Initial::Constract(Eigen::MatrixXd& s){
	if(type=="type1"){
		set_type1(s, amplitude);
	}else if(type=="type2"){
		for (int i=0;i<defect_num;i++){
			set_type2(s, X[i], Y[i]);
		}
	}
}