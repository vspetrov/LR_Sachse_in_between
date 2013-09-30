#include "LR_lattice.h"
#include <unistd.h>
#include <vector>
#include <iostream>

#define SAVE_RST_FILE 1
#define SHOW_PROGRESS 1

const int Size = 40;
double D1 = 0.002;
double D2 = 0.5;
double D3 = 0.;
double PacePeriod = 500.0; //milliseconds
void Init_system(double **V, double **Vc, LR_vars **LR, Fibroblast **FB, int **type)
{
	*V   = new double[Size];
	*Vc   = new double[Size];
	*LR  = new LR_vars[Size];
	*FB  = new Fibroblast[Size];

	*type  = new int[Size];


    /* Types:
     * 0 - myocyte
     * 1 - FB
     * 2 - extra myo
     */
    for (int i=0; i<Size/2; i++){
        (*type)[i] = (double)rand()/(double)RAND_MAX > 0.5 ? 1 : 0;
        std::cout << (*type)[i] << " " ;
    }
    std::cout << std::endl;

    for (int i=Size/2; i<Size; i++){
        (*type)[i] = 2;
    }

	//(*type)[1] = 1; //fibroblast

	for (int i=0; i<Size; i++)
	{
		(*Vc)[i]  = 0.   	   ;
        if ((*type)[i] == 0 || (*type)[i] == 2)
		{
            (*V)[i]   = -72.	   ;

			(*LR)[i].m = 0.00231609;
			(*LR)[i].h = 0.973114;
			(*LR)[i].j = 0.84991;
			(*LR)[i].d = 0.00434296;
			(*LR)[i].f = 0.880756;
			(*LR)[i].X = 0.018826;
			(*LR)[i].Cai = 0.000445703;

            if (i < Size/2)
                (*LR)[i].G_K1 = 0.05*(double)rand()/(double)RAND_MAX;
            else
                (*LR)[i].G_K1 = 0.6047  ;
            (*LR)[i].Iext = 0.0;
		}
		else
		{
			//G_B = 6.9e-3
			/*(*V)[i] = -60.;
			(*FB)[i].O_shkr = 0.;
			(*FB)[i].C0 = 9.11e-1;
			(*FB)[i].C1 = 8.57e-2;
			(*FB)[i].C2 = 3.02e-3;
			(*FB)[i].C3 = 4.74e-5;
			(*FB)[i].C4 = 2.79e-7;*/

			//G_B = 10*6.9e-3
			(*V)[i] = -24.5;
			(*FB)[i].C0=0.176258;
			(*FB)[i].C1=0.367597;
			(*FB)[i].C2=0.287479;
			(*FB)[i].C3=0.0998994;
			(*FB)[i].C4=0.0129962;
			(*FB)[i].O_shkr = 0.0555375;
		}
	}


	/*(*V)[0] = 10.;
	(*V)[1] = 10.;
	(*V)[2] = 10.;
	(*V)[3] = 10.;
	(*V)[4] = 10.;*/
}



std::vector<double> SolveEquations(double MaxTime, double *V, double *Vc, LR_vars *LR,  Fibroblast *FB, int *type)
{
    int i; //counting variables
    int time;//time itarator
    int MT;
    MT=int(MaxTime/dt);
    double write_interval = 10.0;
    int writeInterval = (int)(write_interval/dt);
    int writeCounter = 0;

    int progressStep = (int)(MaxTime/100.0/dt);
    std::vector<std::pair<double,double> > amplitude(Size);
    for (int i=0; i<Size; i++){
        amplitude[i].first = 1e10; //min
        amplitude[i].second = -1e10; //max
    }

    //integrating on the interval from time to time+dt/2
    for (i=0; i<Size; i++)
        Vc[i]=dt/2.*Coupling(i,V,type);

    for (i=0; i<Size; i++)
        V[i]+=Vc[i];

#if SAVE_RST_FILE > 0
    int fd = open("rst.bin",O_CREAT | O_RDWR, S_IREAD | S_IWRITE);
#endif
    const int pacingPeriod = (int)(PacePeriod/dt);
    int paceCounter = 1;
    const double IextAmplitude = 0;//50.0;
    const int paceDuration = (int)(5.0/dt);
    std::vector<std::vector<double> > spikeMoments;
    spikeMoments.resize(Size/2);
    const double threshold_low = -50;
    const double threshold_high = -20;
    std::vector<int> flag;
    flag.resize(Size/2);
    const int freqCalcOffset=(int)(500.0/dt);

    for (time=0; time<MT; time++)
    {
        if (time > pacingPeriod*paceCounter){

            LR[0].Iext = IextAmplitude;
            if (time > pacingPeriod*paceCounter+paceDuration){
                paceCounter++;
                LR[0].Iext = 0;
            }
        }


        for (i=0; i<Size; i++)
            if (type[i] != 1)
                OdeSolve(i,V,LR);
            else
                OdeSolve_fib(i,V,FB);


        for (i=0; i<Size; i++)
            Vc[i]=dt*Coupling(i,V,type);

        for (i=0; i<Size; i++)
            V[i]+=Vc[i];
#if SAVE_RST_FILE > 0
        if (time >= writeInterval*writeCounter){
            write(fd,V,Size*sizeof(double));
            writeCounter++;
        }
#endif

#if SHOW_PROGRESS > 0
        if (time/progressStep*progressStep == time){
            printf("Progress: %d%%\n",time/progressStep);
        }
#endif
        for (i=0; time > 5000 && i<Size ; i++){
            if (V[i] > amplitude[i].second) amplitude[i].second = V[i];
            if (V[i] < amplitude[i].first ) amplitude[i].first  = V[i];
        }

        if (time > freqCalcOffset){
            for (int j=0; j<Size/2; j++){
                if (0 == flag[j] && V[j] > threshold_high){
                    flag[j] = 1;
                    spikeMoments[j].push_back(time*dt);
                }

                if (1 == flag[j] && V[j] < threshold_low){
                    flag[j]=0;
                }
            }
        }
    }
#if SAVE_RST_FILE > 0
    close(fd);
    printf("writeCounter=%d\n",writeCounter);
#endif


    std::vector<double> freq;
    for (int i=0; i<Size/2; i++){
        std::vector<double> *sm = &spikeMoments[i];
        if (sm->size() >= 2)
            freq.push_back(1000.0/((sm->at(sm->size()-1)-sm->at(0))/(sm->size()-1)));
        else
            freq.push_back(0);
    }
    return freq;
}

void OdeSolve(int i,  double *V, LR_vars *LR)
{
	double vd;
	int kstep=1;
	double delta_t=dt;
	for(int j=0; j<kstep; j++)
	{
		vd=VFunction(i,V,LR);

		if(j==0) //decide on the time substep
		{
			kstep=Substeps(vd);
			delta_t=dt/(double)kstep;
		}
		LR[i].Cai += delta_t*CaiFunction(LR[i].Cai,LR[i].d,LR[i].f,V[i]);
		LR[i].m   = mFunction(V[i],LR[i].m,delta_t);
		LR[i].h   = hFunction(V[i],LR[i].h,delta_t);
		LR[i].j   = jFunction(V[i],LR[i].j,delta_t);
		LR[i].d   = dFunction(V[i],LR[i].d,delta_t);
		LR[i].f   = fFunction(V[i],LR[i].f,delta_t);
		LR[i].X   = XFunction(V[i],LR[i].X,delta_t);
		V[i]   += delta_t*vd;
	}
}

void OdeSolve_fib(int i, double *V, Fibroblast *FB)
{
	double dV, dC0, dC1, dC2, dC3, dC4, dO;
	dV = dt*Vf_function(V[i],FB[i].O_shkr);
	dC0 = dt*C0_function(FB[i].C0,FB[i].C1,V[i]);
	dC1 = dt*C1_function(FB[i].C0,FB[i].C1,FB[i].C2,V[i]);
	dC2 = dt*C2_function(FB[i].C1,FB[i].C2,FB[i].C3,V[i]);
	dC3 = dt*C3_function(FB[i].C2,FB[i].C3,FB[i].C4,V[i]);
	dC4 = dt*C4_function(FB[i].C3,FB[i].C4,FB[i].O_shkr,V[i]);
	dO = dt*O_function(FB[i].C4,FB[i].O_shkr);

	V[i] += dV;
	FB[i].C0 += dC0;
	FB[i].C1 += dC1;
	FB[i].C2 += dC2;
	FB[i].C3 += dC3;
	FB[i].C4 += dC4;
	FB[i].O_shkr += dO;



}

inline int Substeps(double &vd)
{
	const int kmax=100;
	int k;

	const int k0=vd>0. ? 5 : 1;
 	k=k0+(int)fabs(vd);

	return k<kmax ? k : kmax;
}

double Coupling(int i, double *V, int *type)
{
	int ln, rn;
    double C = 0;
    if (type[i] != 2){
        if (i != 0) ln=i-1;
        else ln=i;

        if (i != 3*Size/4-1) rn=i+1;
        else rn=i;

        if (type[i] == 1)
        {
            C += D1*(V[ln]+V[rn]-2*V[i])+
                    D2*(V[i+2*Size/4]-V[i]);
        }
        else if (type[i] == 0)
        {
            C += D1*(V[ln]+V[rn]-2*V[i]);
        }
    }else{
        if (i != 3*Size/4) ln=i-1;
        else ln=i;

        if (i != Size-1) rn=i+1;
        else rn=i;

        C += D3*(V[ln]+V[rn]-2*V[i])+
                D2*(V[i-2*Size/4]-V[i]);
    }
	return C;
}
