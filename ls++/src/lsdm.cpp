#include "libls.hpp"
#ifndef __xlc__
//solve some XLC compiler crazyness...
#include "libls-io.hpp"
#endif

#include <fstream>
#include <time.h>


int main(int argc, char **argv)
{
   if (argc<2) ERROR("Syntax is: lsdm-cheb inputfile" );
#ifdef TB_MPI
    MPI_Init( &argc, &argv );
#endif
    MatrixOptions<CMatrix> cops; MatrixOptions<RMatrix> rops(cops);
    
     /* */
//     enum {
//     lsf_fastpoly=1,
//     lsf_niextrapolation=2,
//     lsf_niwarmstart=4,
//     lsf_autoqbar=8,
//     lsf_optk=16,
//     lsf_autohamspec=32,
//     lsf_niextra2=64,
//     lsf_niextra3=128,
//     lsf_autolength=256, //defined but not used so far....
//     lsf_keepinverse=512
//
//     };
//   Routine for making a complex matrix
//     int sz=5; 
//     CMatrix AH(sz,sz);
//     LSOptions lso;
//     lso.p=10000;
//     lso.beta=10;
//     lso.mu=0;
//     lso.accu=1e-6;
//     lso.flags=lsf_optk | lsf_autoqbar | lsf_niextrapolation ; // the other flags were not supported

//     for (int i=0; i<sz; ++i)
//     for (int j=0; j<sz; ++j)
//     {
//         AH(i,j)=1./(1.+(i-j)*(i-j));
//     }

    //std::ifstream ff("/home/michele/lavoro/michele/ls++/src/input.dat");
    std::ifstream ifile(argv[1]);

    LSHelper LSH;
    IField<LSHelper> ilsh(LSH,"");
    HRTimer hrt;
    ilsh<<ifile;

//     ilsh<<std::cin;
//     std::cerr<<"Finished reading input...\n";

//     LSH.setops(lso);
//     LSH.setham(AH);

//   std::ofstream ofs;
//   ofs.open("guessmat");
//   IField<CMatrix> ifham(value.HAM,"");
//   ofs << ifham;
//   ofs.close();


    hrt.start();
    bool scr = false, sfg = false; //doro
    LSH.compute(scr, sfg); //Doro (see libls.cpp 4 details)
    hrt.stop();
    double dt=hrt;
    std::cerr<<"total time elapsed: "<<dt<<"\n";
    
    //FMatrix<double> fmrho(LSH.RHO);
    //std::cout<<std::string(FMatrix<double>(LSH.RHO));
    std::ofstream os;
    
    std::cerr<<"finally trace is: "<<trace(LSH.RHO)<<"\n";
//    TBBenchmarks.print(std::cerr);
    /*if (LSH.options().rho_file!="")
    {
        os.open(LSH.options().rho_file.c_str());
        if (!os.good()) ERROR("Unable to write rho to file");
        os<<"mode full\n";
        os.precision(12);
        os<<fmrho;
        os.close();
    }*/
    
    return 0;
}
