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
 * #ifdef TB_MPI
 * #define LSMatrix PCrsMatrix<double>
 * #else
 * #define LSMatrix CrsMatrix<double>
 * #endif */

int main(int argc, char **argv)
{ 
#ifdef TB_MPI
    MPI_Init( &argc, &argv );
    int myrank, nprocs; //Doro

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //Doro
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs); //Doro
    MPI_Barrier(MPI_COMM_WORLD);

#endif

/*#include "libls.hpp"
#include "libls-cheb.hpp"
#ifndef __xlc__
#include "libls-io.hpp"
#endif

#include <fstream>
#include <time.h>




int main(int argc, char **argv)
{
#ifdef TB_MPI
    MPI_Init( &argc, &argv );
#endif

// DMT: test addition.

//     RMatrix hs;
//     IField<RMatrix> ifhs(hs,"");
//     std::ifstream ifs;
//     ifs.open("h_spcrs");
//     ifs >> ifhs;

//     FMatrix<double> h(hs);
//     IField<FMatrix<double> > ofh(h,"");
//     std::ofstream ofs;
//     ofs.open("h");
//     ofs << ofh;

//     FMatrix<double> hf;
//     IField<FMatrix<double> > ifh(hf,"");
//     std::ifstream ifs;
//     ifs.open("h");
//     ifs >> ifh;

*/

    std::pout<<"hier sind wir" <<"\n";
    std::ifstream ifs;
std::pout<<"hier sind wir2" <<"\n";
    ifs.open("hplain");
std::pout<<"hier sind wir3" <<"\n";
    int r, c;
    ifs >> r >> c;
    double *h = new double[r*c];
    for (int i = 0; i < r*c; ++i) ifs >> h[i];

    double *rho = new double[r*c];
    MPI_Comm comm = MPI_COMM_WORLD;

    lsdm_cheb_c(1000000, 10.0, -0.988445, 1e-4, 1e-6, r, c, h, rho /*, &comm */);

    FMatrix<double> rf(r, c, rho);
    RMatrix rs(rf, MPI_COMM_WORLD);

    FMatrix<double> hf(r, c, h);
    RMatrix hs(hf, MPI_COMM_WORLD);

    RMatrix rhs(hs);
    mult(rs, hs, rhs);

    std::pout<<"finally trace(RHO) is:      "<<trace(rs)<<"\n";
    std::pout<<"finally trace(RHO.HAM) is: "<<trace(rhs)<<"\n";


/*#ifdef TB_MPI
    MPI_Finalize();
    exit(0);
#endif */
//     end additions


    if (argc<2) ERROR("Syntax is: lsdm-cheb inputfile");

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

    //FMatrix<double> fmrho(LSH.RHO);
    //std::cout<<std::string(FMatrix<double>(LSH.RHO));
    std::ofstream os;

    RMatrix RH(LSH.HAM);
    mult(LSH.RHO,LSH.HAM,RH);
    std::pout<<"finally trace(RHO) is:      "<<trace(LSH.RHO)<<"\n";
    std::pout<<"finally trace(RHO).HAM) is: "<<trace(RH)<<"\n";

/*
    TBBenchmarks.print(std::pout);
    if (LSH.options().rho_file!="")
    {
        os.open(LSH.options().rho_file.c_str());
        if (!os.good()) ERROR("Unable to write rho to file");
        os<<"mode full\n";
        os.precision(12);
        os<<fmrho;
        os.close();
    }
*/
#ifdef TB_MPI
    MPI_Finalize();
    exit(0);
#else
    return 0;
#endif
}
