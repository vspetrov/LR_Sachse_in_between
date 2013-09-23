#include <time.h>
#include "LR_lattice.h"


int main(int argc, char *argv[])
{


	double *V, *Vc;
	LR_vars *LR;
	Fibroblast *FB;




	int *type;


    Init_system(&V,&Vc,&LR,&FB,&type);

    SolveEquations(1000,V,Vc,LR,FB,type);
}
