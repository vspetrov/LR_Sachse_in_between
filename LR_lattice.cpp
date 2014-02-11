#include "LR_lattice.h"
#include <vector>
#include <iostream>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "io.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

#define SAVE_RST_FILE 0
#define SHOW_PROGRESS 1



const int start_myo_size = 0;


double D1 = 0.5;
double D2 = 0;
double D3 = 0.;
double PacePeriod = 500.0; //milliseconds


static double *V;
static int *type;

const int lattice_size=100;
typedef struct obstacle_position_t{
    int i;
    int j;
    int R;
} obstacle_position_t;
const obstacle_position_t obstacle_position = {lattice_size/2,lattice_size/2,
                                               lattice_size/5};

int myo_num=0;
int extra_num=0;

static LR_vars* LR_myo = NULL;
static Fibroblast* FB = NULL;

enum{
    LR_CELL,
    FIB_CELL=255,
    LAST_CELL
};

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
static double *V_draw;
static void displayLattice(){
    memcpy(V_draw,V,sizeof(double)*lattice_size*lattice_size);
    const double Vmin = -80.0;
    const double Vmax = 10.0;
    for (int i=0; i<lattice_size*lattice_size; i++){
        if (V_draw[i] < Vmin) V_draw[i] = Vmin;
        if (V_draw[i] > Vmax) V_draw[i] = Vmax;
    }
    cv::Mat m(cv::Size(lattice_size,lattice_size),CV_64FC1,(void*)V_draw);
    cv::Mat show;
    show.create(m.size(),CV_8UC3);

    for (int i=0; i<show.rows; i++){
        for (int j=0; j<show.cols; j++){
            uchar R,G,B;
            double v = m.at<double>(j,i);
            v = (v-Vmin)/(Vmax-Vmin)*4;
            if (v < 1){
                R = 255;
                B = 0;
                G = (uchar)255*v;
            } else if (v < 2){
                R = (uchar)(255*(2-v));
                B = 0;
                G = 255;
            } else if (v < 3){
                G = 255;
                B = (uchar)255*(v-2);
                R = 0;
            } else {
                G = (uchar)(255*(4-v));
                B = 255;
                R = 0;
            }
            show.at<cv::Vec3b>(j,i)[0] = R;
            show.at<cv::Vec3b>(j,i)[1] = G;
            show.at<cv::Vec3b>(j,i)[2] = B;
        }
    }

    cv:resize(show,show,cv::Size(600,600),cv::INTER_LANCZOS4);
    cv::blur(show,show,cv::Size(10,10));
    cv::imshow("show",show);
    cv::waitKey(1);
}


static void saveState(char *filename){
    int fd = open(filename, O_CREAT | O_RDWR | O_BINARY, S_IREAD | S_IWRITE);
    write(fd,V,sizeof(double)*lattice_size*lattice_size);
    write(fd,LR_myo,sizeof(LR_vars)*myo_num);
    write(fd,FB,sizeof(Fibroblast)*extra_num);
    fd = close(fd);
}

static void readState(char *filename){
    int fd = open(filename, O_RDWR| O_BINARY, S_IREAD | S_IWRITE);
    read(fd,V,sizeof(double)*lattice_size*lattice_size);
    read(fd,LR_myo,sizeof(LR_vars)*myo_num);
    read(fd,FB,sizeof(Fibroblast)*extra_num);
    fd = close(fd);
    for (int i=0; i<myo_num; i++){
        LR_myo[i].V = &V[LR_myo[i].id];
    }
    for (int i=0; i<extra_num; i++){
        FB[i].V = &V[FB[i].id];
        FB[i].extra.V = &FB[i].V_extra;
    }
}


void Init_system()
{
    std::cout << "Entering Init system" << std::endl;
    type = new int[lattice_size*lattice_size];
    for (int i=0; i<lattice_size; i++){
        for (int j=0; j<lattice_size; j++){
            if (((i-obstacle_position.i)*(i-obstacle_position.i)+
                 (j-obstacle_position.j)*(j-obstacle_position.j) <
                 obstacle_position.R*obstacle_position.R)){
                type[i*lattice_size+j] = FIB_CELL;
                extra_num++;
            } else {
                type[i*lattice_size+j] = LR_CELL;
                myo_num++;
            }
        }
    }

    V = new double[lattice_size*lattice_size];
    V_draw = new double[lattice_size*lattice_size];
    FB = new Fibroblast[extra_num];
    LR_myo = new LR_vars[myo_num];
    int myo_counter = 0, extra_counter = 0;
    for (int i=0; i<lattice_size*lattice_size; i++){
        if (type[i] == LR_CELL){
            V[i] = -80;
            setLR(&LR_myo[myo_counter]);
            LR_myo[myo_counter].V = &V[i];
            LR_myo[myo_counter].id = i;
            myo_counter++;
        }else{
            V[i] = -60;
            setFB(&FB[extra_counter]);
            FB[extra_counter].V = &V[i];
            FB[extra_counter].id = i;
            FB[extra_counter].V_extra = -80;
            setLR(&FB[extra_counter].extra);
            FB[extra_counter].extra.V = &FB[extra_counter].V_extra;
            extra_counter++;
        }
    }
    for (int i=0; i<lattice_size; i++){
        V[i] = -30;
    }
}

void cleanUp(){
    delete[] type;
}

struct parallelCoupling
{
    double *dV;
    parallelCoupling(double *DV) : dV(DV) {}

    void operator()(const tbb::blocked_range<size_t>& r) const {
        for (size_t k=r.begin();k!=r.end();++k){
                int i,j;
            int ln, rn, tn, bn;
            i = k / lattice_size;
            j = k % lattice_size;
            ln = i - 1 < 0 ? 0 : i - 1;
            rn = i + 1 > lattice_size - 1? lattice_size - 1 : i + 1;
            tn = j - 1 < 0 ? 0 : j - 1;
            bn = j + 1 > lattice_size - 1? lattice_size - 1 : j + 1;
            dV[i*lattice_size+j] =
                dt*D1*(V[ln*lattice_size+j]+V[rn*lattice_size+j]
                       +V[i*lattice_size+tn]+V[i*lattice_size+bn]
                       - 4*V[i*lattice_size+j]);

        }
    }
};

struct parallelODEsolveMyo
{
    double *dV;
    parallelODEsolveMyo(double *DV) : dV(DV) {}

void operator()(const tbb::blocked_range<size_t>& r) const {
    for (size_t i=r.begin();i!=r.end();++i){
        OdeSolve(LR_myo[i].V,&LR_myo[i]);
        *LR_myo[i].V += dV[LR_myo[i].id];
    }
}
};


struct parallelODEsolveFib
{
    double *dV;
    parallelODEsolveFib(double *DV) : dV(DV) {}

void operator()(const tbb::blocked_range<size_t>& r) const {
    for (size_t i=r.begin();i!=r.end();++i){
        LR_vars *lr_extra = &FB[i].extra;
        double dv = dt*D2*(*lr_extra->V - *FB[i].V);
        OdeSolve_fib(FB[i].V,&FB[i]);
        OdeSolve(lr_extra->V,lr_extra);
        lr_extra->dV = -dv;
        FB[i].dV = dv;

        *FB[i].V += FB[i].dV + dV[FB[i].id];
        *lr_extra->V += lr_extra->dV;

    }
}
};

int SolveEquations(double MaxTime)
{
    int time;//time itarator
    int MT;
    MT=int(MaxTime/dt);

    double write_interval = 2.0;
    int writeInterval = (int)(write_interval/dt);
    int writeCounter = 0;
    int progressStep = (int)(MaxTime/100.0/dt);

    double *dV = new double[lattice_size*lattice_size];

    readState("state.bin");

    tbb::task_scheduler_init init(8);

    for (time=0; time<MT; time++)
    {
        tbb::parallel_for(tbb::blocked_range<size_t>(0,lattice_size*lattice_size),
                          parallelCoupling(dV));

        tbb::parallel_for(tbb::blocked_range<size_t>(0,myo_num),
                          parallelODEsolveMyo(dV));

        tbb::parallel_for(tbb::blocked_range<size_t>(0,extra_num),
                          parallelODEsolveFib(dV));

#if SHOW_PROGRESS > 0
        if (time/progressStep*progressStep == time){
            printf("Progress: %d%%\n",time/progressStep);fflush(stdout);
        }
        displayLattice();
#endif

    }
    // saveState("state.bin");
    cleanUp();
    return 0;
}

void OdeSolve(double *V, LR_vars *LR)
{
    double vd;
    int kstep=1;
    double delta_t=dt;
    for(int j=0; j<kstep; j++)
    {
        vd=VFunction(*V,LR);

        if(j==0) //decide on the time substep
        {
            kstep=Substeps(vd);
            delta_t=dt/(double)kstep;
        }
        LR->Cai += delta_t*CaiFunction(LR->Cai,LR->d,LR->f,*V);
        LR->m   = mFunction(*V,LR->m,delta_t);
        LR->h   = hFunction(*V,LR->h,delta_t);
        LR->j   = jFunction(*V,LR->j,delta_t);
        LR->d   = dFunction(*V,LR->d,delta_t);
        LR->f   = fFunction(*V,LR->f,delta_t);
        LR->X   = XFunction(*V,LR->X,delta_t);
        *V   += delta_t*vd;
    }
}

void OdeSolve_fib(double *V, Fibroblast *FB)
{
    double dV, dC0, dC1, dC2, dC3, dC4, dO;
    dV = dt*Vf_function(*V,FB->O_shkr);
    dC0 = dt*C0_function(FB->C0,FB->C1,*V);
    dC1 = dt*C1_function(FB->C0,FB->C1,FB->C2,*V);
    dC2 = dt*C2_function(FB->C1,FB->C2,FB->C3,*V);
    dC3 = dt*C3_function(FB->C2,FB->C3,FB->C4,*V);
    dC4 = dt*C4_function(FB->C3,FB->C4,FB->O_shkr,*V);
    dO = dt*O_function(FB->C4,FB->O_shkr);

    *V += dV;
    FB->C0 += dC0;
    FB->C1 += dC1;
    FB->C2 += dC2;
    FB->C3 += dC3;
    FB->C4 += dC4;
    FB->O_shkr += dO;



}

inline int Substeps(double &vd)
{
    const int kmax=100;
    int k;

    const int k0=vd>0. ? 5 : 1;
    k=k0+(int)fabs(vd);

    return k<kmax ? k : kmax;
}
