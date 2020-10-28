#ifndef __LIBLS_HPP
#define __LIBLS_HPP 0

/* for some - weird - reason icc complains if we put this after including ioparser*/
typedef enum { 
    lsf_fastpoly=1, 
    lsf_niextrapolation=2, 
    lsf_niwarmstart=4, 
    lsf_autoqbar=8, 
    lsf_optk=16, 
    lsf_autohamspec=32, 
    lsf_niextra2=64,
    lsf_niextra3=128,
    lsf_autolength=256, //defined but not used so far....
    lsf_keepinverse=512

} LSFMode; 

#include <complex> 
#include "tbdefs.hpp"
#include "matrix-crs.hpp"
#include "matrix-pcrs.hpp"
#include "matrix-full.hpp"
#include "matrix-full-blas.hpp"
#include "matrix-fun.hpp" 
#include "minsearch.hpp"
#include "matrix-conv.hpp"
#include "matrix-io.hpp"
#include "ioparser.hpp"
 
namespace toolbox {
     std::istream& operator>>(std::istream& istr, LSFMode& lsf);
     std::ostream& operator<<(std::ostream& ostr, const LSFMode& lsf);
 };

using namespace toolbox;

typedef std::complex<double> stdcomplex;
#ifdef TB_MPI
#define Matrix CrsMatrix
#define Matrix2 FMatrix
#else
#define Matrix CrsMatrix
#endif

typedef Matrix<stdcomplex> CMatrix;
typedef Matrix2<stdcomplex> C2Matrix; 
typedef Matrix<double> RMatrix;
  


class LSHelper;
typedef struct _LSOptions{
    unsigned long p;
    unsigned long qbar;
    unsigned long flags;
    double beta, mu;             //1/kbT, chem.potential 
    double hs_max, hs_min;        //boundaries of ham. spectrum
    double accu;                 //accuracy to be seeked for
    bool frestart;
    std::string ham_file;
    std::string rho_file; 
    MatrixOptions<CMatrix> mops;
    _LSOptions(): p(1), qbar(1), flags(0), beta(1.), mu(0.), hs_max(0.), hs_min(0.), accu(1e-4), frestart(false), mops() {}
} LSOptions;

//function to compute convergency ratio for series expansion:
//the two parameters are r and th in k=1+r/p exp[i th] form
class fscratio {
    private: 
        std::valarray<double> x;
//!this must be remade private!
        public:        double q, p, hs_max, hs_min; //fixed par
        double vq, wq, sp, sm, p2;
    public:
        void set_pars(double np, double nq, double nhs_max, double nhs_min)
        {
            q=nq; p=np; hs_max=nhs_max; hs_min=nhs_min;
            double t1;
            t1=constant::pi *(q+q-1)/(p+p);
            vq=cos(t1); wq=sin(t1);
            sp=exp(-hs_max/(p+p)); sm=exp(-hs_min/(p+p));
            sp-=1; sm-=1;
            p2=p*p;
        }
        
        unsigned long size() { return 2; } 
        void set_vars(const std::valarray<double>& vars) { x=vars; }
        void get_vars(std::valarray<double>& vars) const {if (vars.size()!=2) vars.resize(2); vars=x;}
        void get_value(double &rv)
        {
            double r2=x[0]*x[0],pr=p*x[0], ct=cos(x[1]), st=sin(x[1]);
        
            double den,rvp,rvm;
            den=r2-2*p2*(vq-1)+2* pr* (wq*st-(vq-1)*ct);
            rvp=(r2+p2*sp*sp-2*pr*sp*ct)/den;
            rvm=(r2+p2*sm*sm-2*pr*sm*ct)/den;
            rv=(rvp>rvm?rvp:rvm);
        }

        bool bailout() {return false;}
        fscratio() : x((double) 0.,2) {}
};

class LSHelper
{
    friend class toolbox::IField<LSHelper>;
    friend class toolbox::IField<const LSHelper>;
private: 
    LSOptions pops; 
    CMatrix HAM;
    RMatrix ETHAMP;
    bool finit;
    void init_all(); 
    public:
    RMatrix RHO;
    void opt_k(const double& q, stdcomplex& k, double &csi);
    void opt_qbar();
    void seriescoeff(const double &q, const stdcomplex& k, std::valarray<stdcomplex>& sc, const double& scale=1.);
    void tailcoeff(const double &q, const stdcomplex& k, std::valarray<stdcomplex>& tc, const double& scale=1.);

public:
    LSHelper() : finit(false) {};
    LSHelper(const LSOptions& nops) : finit(false) {setops(nops);}
    const LSOptions& options() {return pops;}
    void setops(const LSOptions& nops);
    void setham(const RMatrix& nham);
    void compute(bool& scratch, bool& safeguess); //mod by Doro
};

/*
#ifndef __xlc__
namespace toolbox{
    template<> bool IField<LSOptions>::operator>>(std::ostream& ostr) const; 
    template<> bool IField<const LSOptions>::operator>>(std::ostream& ostr) const; 
    template<> bool IField<LSOptions>::operator<<(std::istream& istr); 
}

namespace toolbox{
    template<> bool IField<LSHelper>::operator>>(std::ostream& ostr) const; 
    template<> bool IField<const LSHelper>::operator>>(std::ostream& ostr) const; 
    template<> bool IField<LSHelper>::operator<<(std::istream& istr); 
}
#else*/
//XLC gets mad about these specialized functions, since he wants them desperately to be inline.
//we make it happy by defining them here.
#ifdef __xlc__
#include "libls-io.hpp"
#endif

#endif //ends ifndef __LIBLS_HPP
