#ifndef _CALCDIF_H_
#define _CALCDIF_H_
#include "ioput.h"

class Calcdif
{
private:
    std::string FilePath;
    
    double t_total;

public:
    double dt;
    int Nt;
    int N_output;
    double dx;
    int Nx;
    double u_outside;

    void readinput();
    double Dx2(int p, int q, Eigen::MatrixXd& s);
    double Dy2(int p, int q, Eigen::MatrixXd& s);
    double Dx4(int p, int q, Eigen::MatrixXd& s);
    double Dy4(int p, int q, Eigen::MatrixXd& s);
    double Dx2y2(int p, int q, Eigen::MatrixXd& s);

    double laplacian(int p, int q, Eigen::MatrixXd& s);
    double laplacian2(int p, int q, Eigen::MatrixXd& s);

    Eigen::MatrixXd grad_x(Eigen::MatrixXd& s);
	Eigen::MatrixXd grad_y(Eigen::MatrixXd& s);
    Eigen::MatrixXd theta(Eigen::MatrixXd& gx, Eigen::MatrixXd& gy);
};
#endif