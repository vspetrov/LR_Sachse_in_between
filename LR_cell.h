#ifndef LR_cell_H
#define LR_cell_H
#include "math.h"
#include <stdio.h>

struct LR_vars{
	double m;
	double h;
	double j;
	double d;
	double f;
	double X;
	double Cai;
	double G_K1;
    double Iext;
};
//right-hand parts of the Luo-Rudy System
extern double VFunction(int i, double *V, LR_vars *LR);
extern double mFunction(double V, double mG, double &delta_t);
extern double hFunction(double V, double hG, double &delta_t);
extern double jFunction(double V, double jG, double &delta_t);
extern double dFunction(double V, double dG, double &delta_t);
extern double fFunction(double V, double fG, double &delta_t);
extern double XFunction(double V, double XG, double &delta_t);
extern double CaiFunction(double Cai,double dG, double fG, double V);
#endif
