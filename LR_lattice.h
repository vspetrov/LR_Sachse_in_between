#ifndef LR_LATTICE_H
#define LR_LATTICE_H
#include "LR_cell.h"
#include "Sachse_fibroblast.h"
#include "stdlib.h"
#include "fcntl.h"
#include "sys/stat.h"

#include <string>



extern void Init_system();//initializing function
void OdeSolve(double *V, LR_vars *LR);
void OdeSolve_fib(double *V, Fibroblast *FB);
inline int Substeps(double &vd);//devides step length due to value of Voltage (V)
int SolveEquations(double MaxTime);




const double dt=0.1;
extern double D1;
extern double D2;
extern double D3;
extern double PacePeriod;
#endif
