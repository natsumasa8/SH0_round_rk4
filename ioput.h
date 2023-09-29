#ifndef _IOPUT_H_
#define _IOPUT_H_

#include <string>

class IOput
{
public:
    int Nx;
    Eigen::MatrixXd input_csv(char* filename);
    void output(Eigen::MatrixXd s, Eigen::MatrixXd s_theta, Eigen::MatrixXd s_energy, char* filename);
    void output_csv(Eigen::MatrixXd s, char* filename);
    void read_vtk(filename);
};

#endif