#ifndef _INITIAL_H_
#define _INITIAL_H_

class Initial
{
	private:
	std::string type; //type1:random, type2:vortex
	int defect_num;
	int X[4], Y[4];
	int r;
	double amplitude;

	public:
	void readinitial();
	void set_type1(Eigen::MatrixXd& s, double amplitude); //random
	void set_type2(Eigen::MatrixXd& s, int X, int Y); //vortex
	void Constract(Eigen::MatrixXd& s);
};
#endif