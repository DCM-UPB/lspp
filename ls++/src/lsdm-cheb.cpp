#include "libls.hpp"
#include "libls-cheb.hpp"
#ifndef __xlc__
#include "libls-io.hpp"
#endif

#include <fstream>
#include <time.h>
#ifdef TB_MPI
#include "mpi.h"
#endif

using namespace toolbox;
#ifdef TB_MPI
#define LSMatrix PCrsMatrix<double>
#else
#define LSMatrix CrsMatrix<double>
#endif


/* using namespace toolbox;
#ifdef TB_MPI
#define LSMatrix PCrsMatrix<double>
#else
#define LSMatrix CrsMatrix<double>
#endif */

int main(int argc, char **argv)
{ 
#ifdef TB_MPI
    MPI_Init( &argc, &argv );
    int myrank, nprocs; //Doro
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //Doro
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //Doro
    MPI_Barrier(MPI_COMM_WORLD);

#endif

    if (argc<2) ERROR("Syntax is: lsdm-cheb inputfile" );
    
    //we replicate input on all nodes, this should not be too bad, as the matrices IO is done better,
    //one node reads then spreads the word
    std::ifstream ifile(argv[1]);
   
    LSChebHelper LSH;

    IField<LSChebHelper> ilsh(LSH,"");
    
    ilsh<<ifile;
    std::perr<<"Finished reading input...\n";
        

    HRTimer hrt;
    hrt.start();
    bool scr = false, sfg = false; //Doro
    LSH.compute(scr, sfg); //Doro - two booleans added - see libls-cheb.cpp 4 details
    hrt.stop();            
    double dt=hrt;
    std::perr<<"total time elapsed: "<<dt<<"\n";
    
    FMatrix<double> fmrho(LSH.RHO);
    //std::cout<<std::string(FMatrix<double>(LSH.RHO));
    std::ofstream os;
    
    RMatrix RH(LSH.HAM);
    mult(LSH.RHO,LSH.HAM,RH);
    std::pout<<"finally trace(RHO) is:      "<<trace(LSH.RHO)<<"\n";
    std::pout<<"finally trace(RHO).HAM) is: "<<trace(RH)<<"\n";
    

    //TBBenchmarks.print(std::pout);
    if (LSH.options().rho_file!="")
    {
        os.open(LSH.options().rho_file.c_str());
        if (!os.good()) ERROR("Unable to write rho to file");
        os<<"mode full\n";
        os.precision(12);
        os<<fmrho;
        os.close();
    }


#ifdef TB_MPI
   
    MPI_Finalize();
    exit(0);
#else
    return 0;
#endif
}
