#ifndef _SHSIMULATE_H_
#define _SHSIMULATE_H_

#include "Calcdif.h"

class SHsimulate : public Calcdif
{
private:
    double a, b, c;
    std::string FilePath;

// protected:
//     Calcdif Cd;

public:
    void ConstractParameters();
    Eigen::MatrixXd CalcSH(Eigen::MatrixXd& s, Calcdif& clc);
    Eigen::MatrixXd CalcRX(Eigen::MatrixXd& s, Calcdif& clc);
    void round_boundary(Eigen::MatrixXd& s, double b);
    Eigen::MatrixXd CalcED(Eigen::MatrixXd& s, Eigen::MatrixXd& gs_x, Eigen::MatrixXd& gs_y, Calcdif& clc); //Energy density
    double CalcE(Eigen::MatrixXd& s, Calcdif& clc); //total energy
    void output_energy(int step, double energy, const char* filename, Calcdif& clc);
};

#endif