#include "LR_lattice.h"
#include <unistd.h>
#include <vector>
#include <iostream>
#define SAVE_RST_FILE 1
#define SHOW_PROGRESS 1



const int start_myo_size = 30;
const int center_fib_size = 100;
const int end_myo_size = 30;
const int extra_myo_size = center_fib_size;

double D1 = 0.8;
double D2 = 0;
double D3 = 0.;
double PacePeriod = 500.0; //milliseconds

static double *V_myo_start = NULL;
static double *V_myo_end = NULL;
static double *V_myo_extra = NULL;
static double *V_fib = NULL;


static LR_vars* LR_myo_start = NULL;
static LR_vars* LR_myo_end = NULL;
static LR_vars* LR_myo_extra = NULL;
static Fibroblast* FB = NULL;

static void setLR(LR_vars *lr){
    lr->m = 0.00231609;
    lr->h = 0.973114;
    lr->j = 0.84991;
    lr->d = 0.00434296;
    lr->f = 0.880756;
    lr->X = 0.018826;
    lr->Cai = 0.000445703;
    lr->G_K1 = 0.6047;
    lr->Iext = 0.0;
}

static void setFB(Fibroblast *fb){
    //G_B = 6.9e-3
    fb->O_shkr = 0.;
    fb->C0 = 9.11e-1;
    fb->C1 = 8.57e-2;
    fb->C2 = 3.02e-3;
    fb->C3 = 4.74e-5;
    fb->C4 = 2.79e-7;

    //G_B = 10*6.9e-3
//    fb->C0=0.176258;
//    fb->C1=0.367597;
//    fb->C2=0.287479;
//    fb->C3=0.0998994;
//    fb->C4=0.0129962;
//    fb->O_shkr = 0.0555375;

    fb->Iext = 0;
}

static void initMyo(double *V, LR_vars *lr, int size){
    for (int i=0; i<size; i++){
        V[i] = -80;
        setLR(&lr[i]);
    }
}
static void initFib(double *V, Fibroblast *fb, int size){
    for (int i=0; i<size; i++){
        V[i] = -60;
        setFB(&fb[i]);
    }
}
void Init_system()
{
    V_myo_start   = new double[start_myo_size];
    V_myo_end   = new double[end_myo_size];
    V_myo_extra   = new double[extra_myo_size];
    V_fib   = new double[center_fib_size];

    LR_myo_start  = new LR_vars[start_myo_size];
    LR_myo_end  = new LR_vars[end_myo_size];
    LR_myo_extra  = new LR_vars[extra_myo_size];

    FB  = new Fibroblast[center_fib_size];

    initMyo(V_myo_start,LR_myo_start,start_myo_size);
    initMyo(V_myo_end,LR_myo_end,end_myo_size);
    initMyo(V_myo_extra,LR_myo_extra,extra_myo_size);

    initFib(V_fib,FB,center_fib_size);
}

void cleanUp(){
    delete[] V_myo_start;
    delete[] V_myo_end;
    delete[] V_myo_extra;
    delete[] V_fib;

    delete[] LR_myo_start;
    delete[] LR_myo_end;
    delete[] LR_myo_extra;
    delete[] FB;
}


std::pair<double, double> SolveEquations(double MaxTime)
{
    int time;//time itarator
    int MT;
    MT=int(MaxTime/dt);
    double write_interval = 2.0;
    int writeInterval = (int)(write_interval/dt);
    int writeCounter = 0;

    int progressStep = (int)(MaxTime/100.0/dt);

    double *dV_myo_start = new double[start_myo_size];
    double *dV_myo_end   = new double[end_myo_size];
    double *dV_myo_extra   = new double[extra_myo_size];
    double *dV_fib   = new double[center_fib_size];


#if SAVE_RST_FILE > 0
    int fd = open("rst.bin",O_CREAT | O_RDWR, S_IREAD | S_IWRITE);
#endif
    const int pacingPeriod = (int)(PacePeriod/dt);
    int paceCounter = 1;
    const double IextAmplitude = 50.0;
    const int paceDuration = (int)(5.0/dt);

    const double threshold_low = -50;
    const double threshold_high = -20;
    int flag = 0;
    const int freqCalcOffset=(int)(500.0/dt);
    std::vector<double> spikeMoments;
    std::pair<double,double> Amplitude(-1e10,1e10);

//    FILE *ofs = fopen("ts_m10.txt","w");
    for (time=0; time<MT; time++)
    {
        if (time > pacingPeriod*paceCounter){

//            FB[0].Iext = IextAmplitude;
            LR_myo_start[0].Iext = IextAmplitude;
            if (time > pacingPeriod*paceCounter+paceDuration){
                paceCounter++;
//                FB[0].Iext = 0;
                LR_myo_start[0].Iext = 0;
            }
        }

        if (time == MT/2){
            D2 = 1.5;
        }

        if (time > (100/dt)){
//            fprintf(ofs,"%g %g\n",time*dt-100,V_fib[0]);
            double v = V_fib[0];
            if (v > Amplitude.first){
                Amplitude.first = v;
            }
            if (v < Amplitude.second){
                Amplitude.second = v;
            }
        }


        for (int i=0; i<start_myo_size; i++){
            OdeSolve(i,V_myo_start,LR_myo_start);
        }
        for (int i=0; i<end_myo_size; i++){
            OdeSolve(i,V_myo_end,LR_myo_end);
        }
        for (int i=0; i<extra_myo_size; i++){
            OdeSolve(i,V_myo_extra,LR_myo_extra);
        }
        for (int i=0; i<center_fib_size; i++){
            OdeSolve_fib(i,V_fib,FB);
        }

        for (int i=0; i<start_myo_size; i++){
            dV_myo_start[i] = 0;
            if (i > 0) dV_myo_start[i] += D1*(V_myo_start[i-1]-V_myo_start[i]);
            if (i <start_myo_size-1) dV_myo_start[i] += D1*(V_myo_start[i+1]-V_myo_start[i]);
            if (i == start_myo_size-1) dV_myo_start[i] += D1*(V_fib[0]-V_myo_start[i]);
        }
        for (int i=0; i<end_myo_size; i++){
            dV_myo_end[i] = 0;
            if (i > 0) dV_myo_end[i] += D1*(V_myo_end[i-1]-V_myo_end[i]);
            if (i <end_myo_size-1) dV_myo_end[i] += D1*(V_myo_end[i+1]-V_myo_end[i]);
            if (i == 0) dV_myo_end[i] += D1*(V_fib[center_fib_size-1]-V_myo_end[i]);
        }
        for (int i=0; i<extra_myo_size; i++){
            dV_myo_extra[i] = D2*(V_fib[i]-V_myo_extra[i]);
            if (i > 0) dV_myo_extra[i] += D3*(V_myo_extra[i-1]-V_myo_extra[i]);
            if (i <extra_myo_size-1) dV_myo_extra[i] += D3*(V_myo_extra[i+1]-V_myo_extra[i]);
        }

        for (int i=0; i<center_fib_size; i++){
            dV_fib[i] = D2*(V_myo_extra[i] - V_fib[i]);
            if (i > 0) dV_fib[i] += D1*(V_fib[i-1]-V_fib[i]);
            if (i <center_fib_size-1) dV_fib[i] += D1*(V_fib[i+1]-V_fib[i]);
            if (start_myo_size > 0)
                if (i == 0) dV_fib[i] += D1*(V_myo_start[start_myo_size-1]-V_fib[i]);
            if (end_myo_size > 0)
                if (i == center_fib_size-1) dV_fib[i] += D1*(V_myo_end[0]-V_fib[i]);
        }

        for (int i=0; i<start_myo_size; i++){
            V_myo_start[i] += dt*dV_myo_start[i];
        }
        for (int i=0; i<end_myo_size; i++){
            V_myo_end[i] += dt*dV_myo_end[i];
        }
        for (int i=0; i<extra_myo_size; i++){
            V_myo_extra[i] += dt*dV_myo_extra[i];
        }
        for (int i=0; i<center_fib_size; i++){
            V_fib[i] += dt*dV_fib[i];
        }


#if SAVE_RST_FILE > 0
        if (time >= writeInterval*writeCounter){
            write(fd,V_myo_start,start_myo_size*sizeof(double));
            write(fd,V_fib,center_fib_size*sizeof(double));
            write(fd,V_myo_end,end_myo_size*sizeof(double));
            writeCounter++;
        }
#endif

#if SHOW_PROGRESS > 0
        if (time/progressStep*progressStep == time){
            printf("Progress: %d%%\n",time/progressStep);
        }
#endif

        if (time > freqCalcOffset){
            double v = V_myo_end[end_myo_size-1];
            if (0 == flag && v > threshold_high){
                flag = 1;
                spikeMoments.push_back(time*dt);
            }

            if (1 == flag && v < threshold_low){
                flag=0;
            }
        }
    }
#if SAVE_RST_FILE > 0
    close(fd);
    printf("writeCounter=%d\n",writeCounter);
#endif

//    fclose(ofs);
    double freq;
    const double refFreq = 1000.0/PacePeriod;
    if (spikeMoments.size() >= 2)
        freq = 1000.0/((spikeMoments[spikeMoments.size()-1]-spikeMoments[0])/(spikeMoments.size()-1));
    else
        freq = 0;

    cleanUp();
    return Amplitude;
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
    dV = dt*Vf_function(V[i],FB[i].O_shkr,FB[i].Iext);
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


