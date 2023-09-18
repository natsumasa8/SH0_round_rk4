#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include "ioput.h"
#include "Calcdif.h"

Eigen::MatrixXd IOput::input_csv(char* filename){
	std::ifstream ifs(filename);
	if (!ifs){
		std::cout << "ioput.cpp:input_csv error: unabel to open " << filename << "." << std::endl;
	}
	Eigen::MatrixXd laodMatrix;
	std::string line;
	int rows = 0;

	while (std::getline(ifs, line)){
		std::stringstream ss(line);
		double value;
		int cols = 0;

		while (ss >> value){
			if (cols == 0){
				laodMatrix.conservativeResize(rows + 1, line.size());
			}
			laodMatrix(rows, cols) = value;
			++cols;

			if (ss.peek() == ','){
				ss.ignore();
			}
		}
		++rows;
	}
	ifs.close();

	std::cout << "CSV " << filename << " was read." << std::endl;
	return laodMatrix;
}


void IOput::output(Eigen::MatrixXd s, Eigen::MatrixXd s_theta, Eigen::MatrixXd s_energy, char* filename){
    std::ofstream ofs(filename);
    if(!ofs){
		std::cout << "ioput.cpp:output error: unable to open " << filename << "." << std::endl;
	}
	ofs << "# vtk DataFile Version 3.0" << std::endl;
	ofs << "SH2vortex" << std::endl;
	ofs << "ASCII" << std::endl;
	ofs << "DATASET STRUCTURED_POINTS" << std::endl;
	ofs << "DIMENSIONS " << Nx << " " << Nx << " " << 1 << std::endl; //Nx * Nx
	ofs << "ORIGIN " << 0 << " " << 0 << " " << 0 << std::endl;
	ofs << "SPACING " << 1 << " " << 1 << " " << 1 << std::endl;
	ofs << "POINT_DATA " << Nx*Nx*1 << std::endl; //Nx * Nx

//vector u size :  Nx * Nx
	ofs << "SCALARS u double" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	ofs << s << std::endl;
//vector u_theta size :Nx * Nx
    ofs << "SCALARS u_theta double" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
	ofs << s_theta << std::endl;
// //vector u_grad size :Nx * Nx
    ofs << "SCALARS u_energy_density double" << std::endl;
	ofs << "LOOKUP_TABLE default" << std::endl;
    ofs << s_energy << std::endl;

	ofs.close();
}

void IOput::output_csv(Eigen::MatrixXd s, char* filename){
	std::ofstream ofs(filename);
	if(!ofs){
		std::cout << "ioput.cpp:output_csv error: unable to open " << filename << "." << std::endl;
	}
	for (int i=0; i<Nx; i++){
		for(int j=0; j<Nx; j++){
			if (j<Nx-1){
				ofs << s(i,j) << ",";
			}else{
				ofs << s(i,j) << "\n";
			}
		}
	}
	std::cout << "csv: " << filename << std::endl;
	ofs.close();
}