#ifndef LR_LATTICE_H
#define LR_LATTICE_H
#include "LR_cell.h"
#include "Sachse_fibroblast.h"
#include "stdlib.h"
#include "fcntl.h"
#include "sys/stat.h"




double Coupling(int i,  double *V, int *type); // This function returns the value of coupling D*(...) of the cell [i,j]
extern void Init_system(double **V, double **Vc, LR_vars **LR, Fibroblast **FB, int **type);//initializing function
void OdeSolve(int i, double *V, LR_vars *LR);
void OdeSolve_fib(int i, double *V, Fibroblast *FB);
inline int Substeps(double &vd);//devides step length due to value of Voltage (V)
extern double SolveEquations(double MaxTime, double *V, double *Vc, LR_vars *LR, Fibroblast *FB, int *type);//Solves the task



extern const int Size;

const double dt=0.1;
extern double D1;
extern double D2;
extern double D3;
extern double PacePeriod;
#endif
