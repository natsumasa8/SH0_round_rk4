#include <iostream>
#include <random>
#include <Eigen/Core>
#include "initial.h"

// void Initial::readinitial(){
// 	std::cout << "Put initial condition file below." << std::endl;

// 	std::string FilePath;
// 	std::cin >> FilePath;
// 	std::ifstream fin(FilePath);
// 	if(!fin){
// 		std::cout << "Error: cannot open input file." << std::endl;
// 		exit(1);
// 	}

// //reading inputfile
// 	std::string dummy, check;
// 	fin >> dummy; fin >> type; std::cout << "type(initial condition): " << type << std::endl;
// 	if(type=="type1"){
// 		//no processing
// 	}elseif(type=="type2"){
// 		fin >> dummy; fin >> r; std::cout << "raduis: " << r << std::endl;
// 		fin >> dummy; fin >> dr; std::cout << "rd: " << dr << std::endl;
// 		fin >> dummy; fin >> N_theta; std::cout << "N_theta" << N_theta << std::endl;
// 	}else{
// 		std::cout << "Error: invalid input of 'type'." << std::endl;
// 	}

// 	fin >> check;
// 	if(check != "End"){
// 		std::cout << "Error: invalid input file." << std::endl;
// 		exit(1);
// 	}
// 	fin.close();
// 	std::cout << "Input has successfuly done." << std::endl;
// }

//random-------------------------------------------------------
void Initial::set_type1(Eigen::MatrixXd& s, double amplitude){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> rand001(-1.0,1.0);

#pragma omp parallel for
    for(int i=0;i<s.rows();i++){
		for(int j=0;j<s.rows();j++){
			s(i,j) += amplitude * rand001(mt);
		}
	}
}
//vortex-------------------------------------------------------
void Initial::set_type2(Eigen::MatrixXd& s, int X, int Y, int r){
	//note : r is pixels not real distance
#pragma omp parallel for
	for(int i=0;i<s.rows();i++){
		for(int j=0;j<s.rows();j++){
			double distance = std::pow(i-X,2) + std::pow(j-Y,2);
			if (((r-1)*(r-1))<= distance && distance<=((r+1)*(r+1))) {
                s(i,j) = 0.1;
			}
		}
	}
}
//-------------------------------------------------
// void Initial::set_type3(){
// 	return 0;
// }
//boundary-------------------------------------------------------
void Initial::set_boundary(int r, Eigen::MatrixXd& s){
	//note : r is pixels not real distance
//#pragma omp parallel for	
    for(int i=0;i<s.rows();i++){
        for(int j=0;j<s.rows();j++){
            // メッシュ座標から円の中心までの距離を計算
            double distance = std::pow(i - s.rows()/2, 2) + std::pow(j - s.rows()/2, 2);
            if (r*r <= distance) {
                s(i,j) = 0.01;
            }
        }
    }
}
//後で
// void Initial::Constract(){
// 	if(type=='type1'){
// 		set_type1();
// 	}elseif(type=='type2'){
// 		set_type2();
// 	}elseif(type=='type3'){
// 		set_type3();
// 	}
// }