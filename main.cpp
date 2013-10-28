#include <time.h>
#include "LR_lattice.h"
#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <sstream>
#include <vector>
#define MPI_2D_SEARCH 0
int main(int argc, char *argv[])
{

#if MPI_2D_SEARCH == 0

    std::vector< std::vector<std::pair<double, double> > > rst;
    std::vector<double> d2;
    for (D2 = 1.35; D2 < 2; D2 += 0.05){
        d2.push_back(D2);
        std::vector<std::pair<double, double> > di_apd_curve;
        for (double delta=0.325; delta <1; delta += 0.025)    {
            Init_system();
            auto DI_APD =  SolveEquations(5000, delta);
            di_apd_curve.push_back(DI_APD);
            //std::cout << "DI: " << DI_APD.first << std::endl;
            //std::cout << "APD: " << DI_APD.second << std::endl;
            //fprintf(ofs,"%g %g\n",DI_APD.first,DI_APD.second);
        }
        rst.push_back(di_apd_curve);
        std::cout << D2 << std::endl;
    }

    FILE *ofs=fopen("di_apd_chain50_d2_varied_pbase500.txt","w");
    fprintf(ofs,"%d ",d2.size());
    for (auto i=0; i<d2.size(); i++)
        fprintf(ofs,"%g ",d2[i]);
    fprintf(ofs,"\n");
    for (auto i=0; i<rst.size(); i++){
        fprintf(ofs,"%d ",rst[i].size());
        for (auto j=0; j<rst[i].size(); j++)
            fprintf(ofs,"%g ", rst[i][j].first);
        fprintf(ofs,"\n");
        fprintf(ofs,"%d ",rst[i].size());
        for (auto j=0; j<rst[i].size(); j++)
            fprintf(ofs,"%g ", rst[i][j].second);
        fprintf(ofs,"\n");
    }
#else
    int size, rank;
    MPI_Init(&argc,&argv);

    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    const double D1min = 0;
    const double D1max = 1.5;
    const double D2min = 0;
    const double D2max = 1.5;
    const int D1steps = 168;
    const int D2steps = 168;
    double *rst = new double[D1steps/size*D2steps];
    for (int i=D1steps/size*rank; i<D1steps/size*(rank+1); i++){
        D1 = D1min + (D1max-D1min)*i/(double)(D1steps-1);
        for (int j=0; j<D2steps; j++){
            D2 = D2min + (D2max-D2min)*j/(double)(D2steps-1);
            double *V, *Vc;
            LR_vars *LR;
            Fibroblast *FB;
            int *type;
            Init_system(&V,&Vc,&LR,&FB,&type);
            double response = SolveEquations(5000,V,Vc,LR,FB,type);
            rst[(i-D1steps/size*rank)*D2steps+j] = response;
            delete[] LR;
            delete[] FB;
            delete[] type;
            delete[] Vc;
            delete[] V;
        }
        if (0==rank){
            std::cout << i+1 << " steps out of " <<D1steps/size << " are done .." << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    double *Allrst = NULL;
    if (rank == 0){
        Allrst = new double[D1steps*D2steps];
    }
    MPI_Gather(rst,D1steps/size*D2steps,MPI_DOUBLE,Allrst,D1steps/size*D2steps,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (rank == 0){
        std::ostringstream oss;
        oss << D3;
        std::string suffix="d3_"+oss.str();
        std::string rst_name = "rst_"+suffix+".bin";
        std::string octave_name = "show_rst_"+suffix+".m";

        int fd = open(rst_name.c_str(),O_CREAT | O_RDWR, S_IREAD | S_IWRITE);
        write(fd,Allrst,D1steps*D2steps*sizeof(double));
        FILE *ofs = fopen(octave_name.c_str(),"w");
        fprintf(ofs,
                "clear all;\n"
                "fd=fopen('%s','r');\n"
                "D1steps=%d;\n"
                "D2steps=%d;\n"
                "D1max=%g;\n"
                "D2max=%g;\n"
                "D1min=%g;\n"
                "D2min=%g;\n"
                "data=fread(fd,[D2steps D1steps],'double');\n"
                "d1=[D1min:(D1max-D1min)/(D1steps-1):D1max];\n"
                "d2=[D2min:(D2max-D2min)/(D2steps-1):D2max];\n"
                "[xx,yy] = meshgrid(d2,d1);\n"
                "surf(xx,yy,data);\n"
                "view(0,90);\n"
                "ylim([D1min D1max]);\n"
                "xlim([D2min D2max]);\n"
                "xlabel('D2');\n"
                "ylabel('D1');\n"
                "shading flat;\n"
                ,rst_name.c_str(),
                D1steps,D2steps,D1max,D2max,D1min,D2min);
        fclose(ofs);
        delete[] Allrst;
        close(fd);
    }
    delete[] rst;
    MPI_Finalize();
#endif
    return 0;
}
