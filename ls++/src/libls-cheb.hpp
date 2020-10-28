#ifndef __LIBLS_CHEB_HPP
#define __LIBLS_CHEB_HPP 0
#include "libls.hpp"
#ifdef TB_MPI
#include <mpi.h>
#endif

using namespace toolbox;

class LSChebHelper
{
    friend class toolbox::IField<LSChebHelper>;
    friend class toolbox::IField<const LSChebHelper>;
    private:


        RMatrix ETHAMP;
        std::valarray<RMatrix> IMq_store;
        bool finit;
        void init_all();
    public:
        LSOptions pops; // DMT: Set public so it is accessible from plain exern "c" function which will be called from fortan.
        RMatrix HAM; //! we make it private or public??? // DMT:
        RMatrix RHO;
        void opt_qbar();
        void chebcoeff(const double &q, const double& z0, const double& zeta, std::valarray<double>& sc);
        void tailchebcoeff(const double &q, const double& z0, const double& zeta, std::valarray<double>& tc);
        unsigned long cheblength(const double& q, const double& z0, const double& zeta);

    public:
        std::ofstream dbgstream; //!DBG
        int myrank;

        LSChebHelper() : finit(false), IMq_store()
        {
#ifdef LSD_DEB
            //!DBG
#ifdef TB_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
            myrank=0;
#endif
            std::string efname;
            efname=std::string("lsdeb.")+int2str(myrank);
            dbgstream.open(efname.c_str(),std::ios_base::out);
            dbgstream<<"DEBUG OUTPUT INITIALIZED\n";
#else
            //debug stream is not initialized, so nothing will be written
#endif
//         MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	std::cerr<<"Node "<<myrank<<" initialized\n";
        };
        ~LSChebHelper() {dbgstream<<"DEBUG CLOSING\n"; dbgstream.close(); }
        LSChebHelper(const LSOptions& nops) : finit(false), IMq_store() {setops(nops);}
        const LSOptions& options() {return pops;}
        void setops(const LSOptions& nops);
        void setham(const RMatrix& nham);
        void compute(bool& scratch, bool& safeguess); //mod by Doro
};
/*
#ifndef __xlc__
namespace toolbox{
    template<> bool IField<LSChebHelper>::operator>>(std::ostream& ostr) const;
    template<> bool IField<const LSChebHelper>::operator>>(std::ostream& ostr) const;
    template<> bool IField<LSChebHelper>::operator<<(std::istream& istr);
}
#else
//XLC gets mad about these specialized functions, since he wants them desperately to be inline.
//we make it happy by defining them here.
#include <libls-io.hpp>
#endif*/
#ifdef __xlc__
#include "libls-io.hpp"
#endif

#endif //ends ifndef __LIBLS_CHEB_HPP

extern "C" void lsdm_cheb_c(int, double, double, double, double, int, int, double*, double* /*, MPI_Comm* */);

