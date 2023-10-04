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
    void ConstructParameters();
    Eigen::MatrixXd SH_rt(Eigen::MatrixXd& s, Calcdif& clc); //Calculating right terms
    Eigen::MatrixXd CalcSH_rk4(Eigen::MatrixXd& s, Calcdif& clc);
    void round_boundary(Eigen::MatrixXd& s, double b);
    Eigen::MatrixXd CalcED(Eigen::MatrixXd& s, Eigen::MatrixXd& gs_x, Eigen::MatrixXd& gs_y, Calcdif& clc); //Energy density
    Eigen::VectorXd CalcED_term(Eigen::MatrixXd& s, Eigen::MatrixXd& gs_x, Eigen::MatrixXd& gs_y, Calcdif& clc);//energy density(term)
    double CalcE(Eigen::MatrixXd& s, Calcdif& clc); //total energy
    void output_energy(int step, double energy, Eigen::VectorXd energy_term, const char* filename, Calcdif& clc);
};

#endif