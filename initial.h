#ifndef _INITIAL_H_
#define _INITIAL_H_

class Initial
{
	private:
	std::string type;

	public:
	void readinitial();
	void set_type1(Eigen::MatrixXd& s); //random
	void set_type2(Eigen::MatrixXd& s, int X, int Y, int r); //vortex
	// void set_type3(); //orientation
	void set_boundary(int r, Eigen::MatrixXd& s);
	// void Constract();
};
#endif