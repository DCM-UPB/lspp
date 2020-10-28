#include "libls-cheb.hpp"
#include <fstream>

void LSChebHelper::setops(const LSOptions& nops)
{
    std::pout<<" * SETTING OPTIONS\n";
    pops=nops;
    HAM.setops(pops.mops);
    if (HAM.rows()>0) init_all();
}

void LSChebHelper::setham(const RMatrix& nham)
{
    std::pout<<" * SETTING THE HAMILTONIAN\n";
    HAM=nham;
    HAM.setops(pops.mops);
    std::pout<<" # HAMILTONIAN DATA: trace: "<<trace(HAM)<<"  norminf: "<<norminf(HAM)<<"  normfrob:"<<normfrob(HAM)<<"\n";
    init_all();
}

void LSChebHelper::init_all()
{
    std::pout<<" * COMPUTING INITIAL MATRICES\n";

    if (pops.flags & lsf_autohamspec)
    {
        std::pout<<" ** EXTIMATING HAMILTONIAN SPECTRUM BOUND\n";
        /*RMatrix tmp,tmp2;
        MatrixOptions<RMatrix> hops; HAM.getops(hops); tmp.setops(hops); tmp2.setops(hops);
        std::<<"now we compute the product\n";

        mult(THAM,THAM,tmp);
        mult(tmp,tmp,tmp2);
       //extimation of the spectral radius of the shifted hamiltonian
        //double x1, x2;
        //x1=srad_ub(tmp2); x2=norminf(tmp2);
        pops.hsr=pow(norminf(tmp2),0.25);

        std::cerr<<"Estimated spectral radius: "<<pops.hsr<<"\n";*/
        /* estimate the bounds by gershgorin's algorithm. we use an auxiliary matrix, as we don't
        want to access directly matrix elements, as we want to work independent of the matrix kind */
        /*
        RMatrix tmp(HAM);
        std::valarray<double> diag(HAM.rows());
        //!probably we should implement a "getdiag" and "setdiag" functions in the matrix library
        for (unsigned long i=0; i<diag.size(); ++i) diag[i]=tmp(i,i);
        tmp-=diag;
        */
        /* dumbest choice: max and min are guessed from a matrix norm */
        pops.hs_min=-(pops.hs_max=norminf(HAM));
    }

    std::pout<<" ** EXTIMATING BEST Q\n";
    if (pops.flags & lsf_autoqbar) opt_qbar();
    std::pout <<"EXTIMATED Q = " << pops.qbar << "\n";
    std::pout<<" ** COMPUTING EXPONENTIALS\n";
    MatrixOptions<RMatrix> mops;
    HAM.getops(mops);
    RMatrix tmp(HAM); tmp-=pops.mu;
    ETHAMP.setops(mops);
    tmp*=(pops.beta/(2*pops.p));
    exp(tmp,ETHAMP, pops.accu*pops.accu);   //computes exponential with improved accuracy: it is cheap, so better being safe than sorry
    std::pout<<"#########################  EXPONENTIAL DATA  #############################\n";
    std::pout<<" trace(1-exp): "<<trace(ETHAMP)-ETHAMP.rows()<<"("<<ETHAMP.rows()<<")"<<"   norminf(1-exp): "<<norminf(ETHAMP)-1.<<" nel: "<<ETHAMP.size()/ETHAMP.rows()<<"\n";

    //for (unsigned long i=0; i<ETHAMP.rows(); ++i)
     //   std::cout <<i<<":"<<ETHAMP(i,i)<<"  ";

    std::pout<<"##########################################################################\n";
    finit=true;
}

//*!! ANALITICAL CHEB COEFFICIENT
void LSChebHelper::chebcoeff(const double &q, const double& z0, const double& zeta, std::valarray<double>& sc)
{
    /*
    std::cerr<<"PARAMETERS: \n"
            <<"q: "<<q<<"\n"
            <<"p: "<<pops.p<<"\n"
            <<"z0: "<<z0<<"\n"
            <<"zeta: "<<zeta<<"\n";
    */
    //computes the cheb. coefficients of N_q^{-1} analitically
    double phiq=(2 * q -1)*constant::pi/(2*pops.p);
    double rz=(cos(phiq)-z0)/zeta,iz=sin(phiq)/zeta;
    stdcomplex w(rz,iz);
    stdcomplex sqw=sqrt(w*w-1.), r, rn;
    if (rz>0) {r=w-sqw; rn=1./sqw;} else {r=w+sqw; rn=-1./sqw;}

//wie: Factor 2 for Chebcoeff 
//    sc=-4./(zeta*zeta*iz);

    sc=-2./(zeta*zeta*iz);
    sc[0]*=0.5*imag(rn);
    for (unsigned long i=1; i<sc.size();  ++i)
    {
        rn*=r;
        sc[i]*=imag(rn);
    }

    /*
    std::ofstream os("ccoeff");
    for (unsigned long i=0; i<sc.size(); ++i) os<<i<<" "<<sc[i]<<"\n";
    os.close();
    */
}
//*/
void LSChebHelper::tailchebcoeff(const double &q, const double& z0, const double& zeta, std::valarray<double>& tc)
{
    //computes the cheb. coefficients of the tail by summing up the analytical contributions.
    std::valarray<double> sc(tc.size());
    //this is already blazing fast, but could be improved by summing up only the relevant terms for large q
    tc=0.;
    for(unsigned long iq=(unsigned long) q; iq<=pops.p; iq++)
    { chebcoeff(iq,z0,zeta,sc); tc+=sc; }
}

#define __CHEB_MAX_LENGTH 100000
unsigned long LSChebHelper::cheblength(const double& q, const double& z0, const double& zeta)
{
    //we want a tight estimate. so we use first a full theoretical estimate, based on a rather
    //rough upper bound, and then we refine it by computing the coefficients for real.
    double xt=ceil(0.5-(pops.hs_max-pops.hs_min)*pops.beta*log(pops.accu)/(constant::pi*(2.*q -1.)));

    //we put a cutoff, in case one puts in unreasonable values
    if (xt> __CHEB_MAX_LENGTH)
    {
        WARNING("Length for Chebyshev series is larger than hardcoded limit.");
        return (unsigned long)  __CHEB_MAX_LENGTH;
    }

    unsigned long lxt=(unsigned long) xt;
    //we compute twice as many terms
    std::valarray<double> tsc(2*lxt+1);
    chebcoeff(q, z0, zeta, tsc);

    //we compute the extra error term (which is going to be tiiiiiny!)
    double phiq=(2 * q -1)*constant::pi/(2*pops.p);
    double rz=(cos(phiq)-z0)/zeta,iz=sin(phiq)/zeta;
    stdcomplex w(rz,iz);
    stdcomplex sqw=sqrt(w*w-1.), r, rn;
    if (rz>0) {r=w-sqw; rn=1./sqw;} else {r=w+sqw; rn=-1./sqw;}
    double err=abs(rn/iz)/(zeta*zeta)*pow(abs(r),2.*xt);

    //now we start summing up terms from the coefficients while we are smaller than the target accuracy
    double target=pops.accu* (4.*pops.p*pops.p)/(constant::pi*constant::pi*(2. *q-1.)*(2.*q-1.));

    unsigned long nt=2*lxt;
    while(err<target && nt > 1) {err+=abs(tsc[nt]); nt--; }

    std::pout<<" * Estimated cheb series length: "<<lxt<<" corrected to: "<<nt<<"\n";
    return nt;
}


//we need a specialized version of newton inverse, because for convergence we better check
//the change in the trace, which is a more reliable estimate of whether the thing has converged
//quindi controlliamo la norma, che dev'essere minore di 1, e poi vediamo che il cambiamento nella
//traccia sia piccolo
bool inverse_newton(
    const RMatrix& M, RMatrix& R,
    IterOptions<double,3> iops=IterOptions<double,3>(__TB_STD_MAXSTEP, __TB_STD_EPS, 0., ichk_default)
)
{
#ifdef BENCHMARK
    unsigned long nm=0;
    TBBenchmarks["LL_inv_newt"].n_calls++; TBBenchmarks["LL_inv_newt"].timer.start();
#endif

    bool fguess=false;
    //we can pass a guess for the inverse in G, if we pass an empty matrix, we start from scratch
    if (R.rows()!=0)
    {
        fguess=true;
    }
    else
    {
        //"guaranteed convergency" initial guess
        std::pout<<"!! guessing safe\n";
        //!transpose(M,R);  this MPI fails... this is not STRICTLY necessary since M is symmetric anyway
        R=M;
        double t=norminf(R)*norminf(M);
        typedef double (*fmap_type) (const double&);
        toolbox::map(R,fmap_type(conj));
        scale(R,1./t);
    }
    RMatrix::index_type n=M.rows();

    MatrixOptions<RMatrix> mopts;
    R.getops(mopts);
    RMatrix tmp=R; tmp.resize(0,0); RMatrix tmp2=tmp;
    //mopts.atthresh*=0.1;  //smaller threshold for intermediate multiply!
    //tmp.setops(mopts);
    double errn;

    std::pout<<" ######################### MATRIX INVERSION #########################\n";
    mult(R,M,tmp);
    std::pout<<" * guess-> trace: "<<trace(R) <<"   ninf: "<<norminf(R)<< " nel: " << R.size()/R.rows()<<std::endl;
    std::pout<<" * emmeq-> trace: "<<trace(M) <<"   ninf: "<<norminf(M)<< " nel: " << M.size()/M.rows()<<std::endl;
    std::pout<<" * first-> trace: "<<trace(tmp) <<"   ninf: "<<norminf(tmp)<< " nel: " << tmp.size()/tmp.rows()<<std::endl;

    while (!iops)
    {
        std::pout<<" * >>> ITERATION ( "<<iops.iter()<<" )  <<< "<<std::endl;

#ifdef BENCHMARK
        nm+=2;
#endif
        mult(tmp,R,tmp2);
        std::pout<<" * RMR -> trace: "<<trace(tmp2) <<"   ninf: "<<norminf(tmp2)<< " nel: " << tmp2.size()/tmp2.rows()<<std::endl;
        scale(R,-2.);
        incr(R, tmp2);
        neg(R);  //here we have the new guess
        std::pout<<" * 2R-RMR -> trace: "<<trace(R) <<"   ninf: "<<norminf(R)<< " nel: " << R.size()/R.rows()<<std::endl;

        mult(R,M,tmp);  //by computing this here, we get the error 'almost' for free
        std::pout<<" * MR  -> trace: "<<trace(tmp) <<"   ninf: "<<norminf(tmp)<< " nel: " << tmp.size()/tmp.rows()<<std::endl;

        tmp-=1.;  //for (typename MC::index_type i=0; i<n; ++i) tmp(i,i)-=1.;
        errn=norminf(tmp);
        //check for divergence
        if (abs(trace(tmp))/n >1.)
        {
            if (!fguess) {ERROR("Matrix inversion is diverging even with 'safe guess'. M should be VERY ill-conditioned");}
            else return false;
            //it is not nice nor useful to fallback without notice to the caller!
/*            R.resize(0,0);
#ifdef BENCHMARK
            TBBenchmarks["LL_inv_newt"].timer.stop(); TBBenchmarks["LL_inv_newt"].tot_time+=TBBenchmarks["LL_inv_newt"].timer;
            TBBenchmarks["LL_inv_newt"].tot_prop["mm_mult"]+=nm;  //n. of multiplications
#endif
            inverse_newton(M,R,iops); //try to invert from scratch
            return;*/
        }
        for (unsigned long i=0; i<2; ++i)
            iops.setval(errn,i);
        iops.setval(trace(R),2);

        //recover the shift
        tmp+=1.; //for (typename MC::index_type i=0; i<n; ++i) tmp(i,i)+=1.;
    }
    std::pout<<" ####################################################################\n";
    if (iops.maxstep_reached()) ERROR("Newton inversion reached max number of steps without converging within desired accuracy\n");
#ifdef BENCHMARK
    TBBenchmarks["LL_inv_newt"].timer.stop(); TBBenchmarks["LL_inv_newt"].tot_time+=TBBenchmarks["LL_inv_newt"].timer;
    TBBenchmarks["LL_inv_newt"].tot_prop["mm_mult"]+=nm;  //n. of multiplications
#endif
    return true;
}


void LSChebHelper::compute(bool& scratch, bool& safeguess)
{
    static int ncalls=0; //!DBG
    if (!finit) init_all();
#ifdef BENCHMARK
    TBBenchmarks["LL_OVERALL"].n_calls++; TBBenchmarks["LL_OVERALL"].timer.start();
#endif
    dbgstream<<ncalls<<" COMPUTE CALL"<<std::endl; ncalls++; //!DBG
    // ***  Declares and set options for working matrices
    RMatrix Nq, INq;
    RHO.resize(HAM.rows(),HAM.rows());
    Nq.resize(HAM.rows(),HAM.rows());
    INq.resize(HAM.rows(),HAM.rows());
    MatrixOptions<RMatrix> rmops;
    HAM.getops(rmops);
    RHO.setops(rmops);
    INq.setops(rmops);
    Nq.setops(rmops);

    double z0, zeta;
    z0=(exp((pops.hs_max-pops.mu)*pops.beta/(2*pops.p))+exp((pops.hs_min-pops.mu)*pops.beta/(2*pops.p)))*0.5;
    //we slightly enlarge zeta to be on the safe side
    zeta=(exp((pops.hs_max-pops.mu)*pops.beta/(2*pops.p))-exp((pops.hs_min-pops.mu)*pops.beta/(2*pops.p)))*0.5*1.01;
    std::pout<<" * Selected values of z0: "<<z0<<" and zeta: "<<zeta<<"\n";

    //COMPUTE CHEBYSHEV LENGTH!!!
    unsigned long ms=cheblength(pops.qbar,z0,zeta);
    std::valarray<double> sc(ms), tc(ms);

    // ******************** does the series part **********************  //
    bool ffastpoly=(pops.flags&lsf_fastpoly), fwarmstart, fhaveprev;
    if (IMq_store.size()==0)
    {
        fhaveprev=false;
        if (pops.flags & lsf_keepinverse) IMq_store.resize(pops.qbar-1);
    }
    else fhaveprev=true;
    fwarmstart=(pops.qbar>1 && (pops.flags & lsf_niwarmstart) && !fhaveprev);
    if (ffastpoly && (fwarmstart?3:2)*sqrt(1. * ms) >= ms)
    {
        ffastpoly=false;
        std::pout<<"Fallback on standard polynomial evaluation, which should require less operations\n";
    }



    Nq=ETHAMP; Nq-=z0; Nq*=1./zeta;
    if (fwarmstart)
    {
        chebcoeff(pops.qbar,z0,zeta,sc);
        tailchebcoeff(pops.qbar,z0,zeta,tc);

        std::valarray<RMatrix *> cm(2); std::valarray<std::valarray<double> *> ck(2);
        cm[0]=&INq;   cm[1]=&RHO;
        ck[0]=&sc;    ck[1]=&tc;

        if (ffastpoly) chebyshev_fastpoly(Nq,ck,cm);
        else chebyshev_poly(Nq,ck,cm);
    }
    else
    {
        tailchebcoeff(pops.qbar,z0,zeta,tc);
        if (ffastpoly) chebyshev_fastpoly(Nq,tc,RHO);
        else chebyshev_poly(Nq,tc,RHO);
        INq.resize(0,0);
    }

    std::pout<<" ###################  CHEB PART DIAGNOSTICS #####################\n";
    std::pout<<" * RHO TRACE: " <<trace(RHO)<<" \n";
    std::pout<<" * RHO NINF:  " <<norminf(RHO)<<" \n";
    std::pout<<" * RHO NFROB: " <<normfrob(RHO)<<" \n";
    if (fwarmstart)
    {
        std::pout<<" * INQ TRACE: " <<trace(INq)<<" \n";
        std::pout<<" * INQ NINF:  " <<norminf(INq)<<" \n";
        std::pout<<" * INQ NFROB: " <<normfrob(INq)<<" \n";
    }
    std::pout<<" ################################################################\n";

    dbgstream<<"CHEBYSHEV DONE!"<<std::endl; //!DBG

    //to the newton part
    /*
    newton inversion bailout options
      we consider inversion converged checking the norm of the residual:
    this NEEDs to be below 1
    if it is below accuracy OR if the change in the norm is below it, we quit
    */

    double nwaccu;
    IterOptions<double,3> inwops(100,pops.accu, 0., ichk_default);
    inwops.thresh[0]=10.; //!MINIMUM we want the residual to be below 1!  DO WE REALLY NEED THIS COND???
    inwops.flags[1]=ichk_sufficient;
    inwops.flags[2]=ichk_change  | ichk_relative;


    std::pout<<" *********************************************************\n";
    std::pout<<" *         COMPUTING NEWTON CONTRIBUTION                 *\n";
    std::pout<<" *********************************************************\n";
    double cq=cos(constant::pi/pops.p), oldcq;
    double kx, kx1, kx2, kx3;
    RMatrix INq_2, INq_3, INq_s; INq_2.setops(rmops); INq_2.resize(0,0); INq_3.setops(rmops); INq_3.resize(0,0);
    //we must compute Nq with hiiigh accuracy, at least as H^2/P^2
    rmops.atthresh*=(pops.hs_max-pops.hs_min)*pops.beta/(pops.p*2.);
    rmops.atthresh*=(pops.hs_max-pops.hs_min)*pops.beta/(pops.p*2.);
    Nq.setops(rmops);
    mult(ETHAMP,ETHAMP,Nq); Nq+=1;
    cq=oldcq=cos(constant::pi*(pops.qbar+pops.qbar-1)/(pops.p+pops.p));
    rmops.atthresh=0.; Nq.setops(rmops);
    incr(Nq,ETHAMP,-2*cq);


    for (unsigned long q=pops.qbar-1; q>0; --q)
    {
        dbgstream<<"+ NEWTON STEP"<<q<<std::endl; //!DBG
        std::pout<<" >>>>>>>>>  computing q-channel: "<<q<<"  <<<<<<<<<<< \n";
        nwaccu=(pops.hs_max-pops.hs_min)/(constant::pi*(2.*q-1.)); nwaccu*=nwaccu; nwaccu+=1.; nwaccu*=pops.accu;
        std::pout<<" * accuracy threshold: "<<nwaccu<<"\n";

        inwops.thresh[1]=nwaccu;
        cq=cos(constant::pi*(q+q-1)/(pops.p+pops.p));
        incr(Nq,ETHAMP,2*(oldcq-cq)); //get the new Nq the cheapest way possible!!!
        oldcq=cq;

        if (!fhaveprev) {
        if (pops.flags & lsf_niextrapolation && !safeguess) {
            kx=1./(cos(-constant::pi/pops.p)-sin(-constant::pi/pops.p) * tan (constant::pi*((q-1)*2-1)/(pops.p*2)));
            INq_s=INq;
            if (INq_2.rows()>0 && (pops.flags & (lsf_niextra2 | lsf_niextra3)))
            {
                //long history extrapolation
                if (INq_3.rows()>0 && (pops.flags & lsf_niextra3))
                {
                    kx=cos(constant::pi*(2.*q-1.)/(pops.p*2.))*2;kx*=kx*kx;
                    kx*=(cos(constant::pi/(2.*pops.p))-cos(constant::pi*(4.*q+7.)/(pops.p*2.)));
                    kx1=2.*cos(constant::pi*(-2. + q)/pops.p) + 2.*cos(q*constant::pi/pops.p) + cos(q*3.*constant::pi/pops.p) -
                            2*cos(3.*constant::pi*(1. + q)/pops.p) - cos(5.*constant::pi*(1. + q)/pops.p) -
                            2*cos(constant::pi*(2. + 3.*q)/pops.p) - 2*cos(constant::pi*(4. + 3.*q)/pops.p) -
                            cos(constant::pi*(3. + 5.*q)/pops.p) - cos(constant::pi*(4. + 5.*q)/pops.p) +
                            cos(constant::pi*(1 - 3.*q)/pops.p) + 2.*cos(constant::pi*(1. -q)/pops.p) +
                            cos(constant::pi*(1. + 3.*q)/pops.p);
                    kx2=2.*cos(3*constant::pi*(1. + q)/pops.p) + cos(5*constant::pi*(1. + q)/pops.p) -
                            2.*cos(constant::pi*(2. + q)/pops.p) - cos(3*constant::pi*(2. + q)/pops.p) -
                            2.*cos(constant::pi*(3. + q)/pops.p) - 2*cos(constant::pi*(4. + q)/pops.p) +
                            2.*cos(constant::pi*(2. + 3.*q)/pops.p) + 2*cos(constant::pi*(4. + 3.*q)/pops.p) -
                            cos(constant::pi*(5. + 3.*q)/pops.p) - cos(constant::pi*(7. + 3*q)/pops.p) +
                            cos(constant::pi*(6. + 5*q)/pops.p) + cos(constant::pi*(7. + 5.*q)/pops.p);
                    kx3=-3.*cos(constant::pi*(-2. + q)/pops.p) - 3.*cos(3*constant::pi*(1. + q)/pops.p) +
                            3.*cos(constant::pi*(2. + q)/pops.p) + 3*cos(constant::pi*(3. + q)/pops.p) -
                            cos(constant::pi*(7. + q)/pops.p) + cos(constant::pi*(7. + 3*q)/pops.p) +
                            cos(constant::pi*(8. + 3.*q)/pops.p) - cos(constant::pi*(8. + 5*q)/pops.p);
                    kx1/=kx;  kx2/=kx;  kx3/=kx;
                    //std::cerr<<"extra second "<<kx1<<" , "<<kx2<<" , "<<kx3<<" \n";
                    INq*=kx1;
                    incr(INq,INq_2,kx2);
                    incr(INq,INq_3,kx3);
                }
                else
                {
                    kx=cos(constant::pi*(2.*q-1.)/(pops.p*2.))*2.; kx*=kx;
                    kx1=(1.+(sin(q*constant::pi/pops.p)+sin(constant::pi*(3.*q+2.)/pops.p)+sin(constant::pi*(3.*q+1.)/pops.p))
                        /sin(constant::pi*(q+1.)/pops.p))/kx;
                    kx2=- (2.*sin(q*constant::pi/pops.p)+sin(constant::pi*(3.*q+3.)/pops.p)-sin(constant::pi*(3.+q)/pops.p))/
                            (kx*sin(constant::pi*(q+1.)/pops.p));
                    //std::cerr<<"extra first "<<kx1<<" , "<<kx2<<" \n";
                    INq*=kx1;
                    INq_2*=kx2;
                    INq+=INq_2;
                }
            }
            else
            {
                //std::cerr<<"extra zeroth "<<kx<<" \n";
                INq*=kx;
            }
            if (pops.flags & lsf_niextra3) INq_3=INq_2;
            if (pops.flags & lsf_niextra2) INq_2=INq_s;
        }
        else
        {
            INq.resize(0,0);
        }
        }
        else INq=IMq_store[q-1];

        bool inres=inverse_newton(Nq,INq,inwops);
        if (inres==false)
        {
            if (pops.flags &lsf_keepinverse)
            {
                if (safeguess) std::cout << "safeguess true" << std::endl;
                if (fhaveprev && !scratch && !safeguess)
                {
                    IMq_store.resize(0);
                    std::perr<<"Inverse failed to converge with stored inverse. Try to compute from scratch.\n";
                    //try to compute without time extrapolation
#ifdef BENCHMARK
                    TBBenchmarks["LL_OVERALL"].timer.stop(); TBBenchmarks["LL_OVERALL"].tot_time+=TBBenchmarks["LL_OVERALL"].timer;
#endif
                    scratch = true;
                    compute(scratch, safeguess);
                    return;
                }
                else if (!fhaveprev && scratch && !safeguess)
                {
                    IMq_store.resize(0);
                    std::perr<<"Inverse failed to converge with stored inverse. Try to compute with safe guess.\n";
                    //try to compute without time extrapolation
#ifdef BENCHMARK
                    TBBenchmarks["LL_OVERALL"].timer.stop(); TBBenchmarks["LL_OVERALL"].tot_time+=TBBenchmarks["LL_OVERALL"].timer;
#endif
                    safeguess = true;
                    compute(scratch, safeguess);
                    return;
                }
                else
                {
                    ERROR("No way to get inverse to converge. Sorry, bad luck!\n");
                }
            }
        }
        if (pops.flags & lsf_keepinverse) IMq_store[q-1]=INq;
        /*!
        std::string fname;
        std::stringstream ss;
        ss<<"invmat-cb-"<<q;
        ss>>fname;
        std::ofstream os(fname.c_str());
        os<< INq;
        os.close();
        */

        RHO+=INq;
        scratch = false; safeguess = false;
        //std::cerr<<"q: "<<q<<" tr(m_q^-1): "<<trace(INq)/(1.*pops.p)<<"\n";
    }

    //reuse the matrices to compute the complete DM
    //!don't mess around: temporary until we implement cheb part
    mult(ETHAMP,ETHAMP,Nq); Nq-=1.;
    mult(Nq,RHO,INq);


    //INq+=(pops.qbar-1.);
    INq*=-0.5/pops.p; INq+=0.5;
    RHO=INq;

//     std::cout << "Tr[rho] : " <<  trace(RHO) << std::endl;
//     std::cout << "Tr[H] : " <<  trace(HAM) << std::endl;

    //stitches the pieces together
    /*RMatrix rRHO; copymap(cRHO,rRHO,re);
    RHO+=rRHO;
    RHO*=-(1./pops.p);
    RHO+=1.;*/
#ifdef BENCHMARK
    TBBenchmarks["LL_OVERALL"].timer.stop(); TBBenchmarks["LL_OVERALL"].tot_time+=TBBenchmarks["LL_OVERALL"].timer;
#endif

}

void LSChebHelper::opt_qbar()
{
    if (pops.flags & lsf_niwarmstart) pops.qbar=2; else pops.qbar=1;
    double z0, zeta;
    z0=(exp((pops.hs_max-pops.mu)*pops.beta/(2*pops.p))+exp((pops.hs_min-pops.mu)*pops.beta/(2*pops.p)))*0.5;
    zeta=(exp((pops.hs_max-pops.mu)*pops.beta/(2*pops.p))-exp((pops.hs_min-pops.mu)*pops.beta/(2*pops.p)))*0.5*1.01;

    double sn, sn_next, nn;
    double deltaseries;
    sn_next=cheblength(pops.qbar,z0,zeta);
    while (pops.qbar<pops.p)
    {
        //compute conv. speed for chebychev
        sn=sn_next; sn_next=cheblength(pops.qbar+1,z0,zeta);

        //compute conv. for newton
        if (pops.flags & lsf_niwarmstart)
        {
            //not so sure....
            nn=2.*ceil(1./log(2.)*log(log(pops.accu)/log(64.*pops.qbar/((1+2.*pops.qbar)*(3+2.*pops.qbar)*(3+2.*pops.qbar)))));
        }
        else
        {
            //boh!
            //cold start extimates
        }
        if (pops.flags & lsf_fastpoly)
        {
            deltaseries=ceil(sqrt(sn))-ceil(sqrt(sn_next));
            if (pops.flags & lsf_niwarmstart) deltaseries*=3.; else deltaseries*=2.;
        }
        else
            deltaseries=ceil(sn)-ceil(sn_next);
        //if (pops.flags & lsf_fastpoly) nm_series=ceil(sqrt(nm_series))*(pops.flags & lsf_niwarmstart?3:2);

        //std::cerr<<"Polynomial order: " << sn<<"\n";
        //std::cerr<<"For qbar= "<<pops.qbar<<" estimated n. newton: "<<nn<<",  n. series: "<<deltaseries<<"\n";
        if (deltaseries<nn) break;
        pops.qbar++;
    }

}



extern "C" void lsdm_cheb_c(int p, double beta, double mu, double accu, double thresh,
                            int r, int c, double *h, double *rho /*, MPI_Comm *comm */) {
// The interface is intentionally basic (i.e. no strucs/classes, templates etc..) to make linking with f08 simple.
// Structs can be made to work too but i'd rather not introduce complexity at the moment (X-2013).
// Much more important issue is the significant amount of copying we do to convert bare c array to pcrs
// It is important to have MPI_Comm* (pointer) as in fortran it is represented by handle. Will be using just mpi_comm_world for now.
//  DMT

    LSOptions o;
    o.p = p;
    o.flags = lsf_fastpoly | lsf_niextrapolation | lsf_niwarmstart | lsf_keepinverse |  lsf_autohamspec | lsf_autoqbar ;
    o.beta = beta;
    o.mu = mu;
    o.accu = accu;
    o.frestart = false;
    o.mops.atthresh = thresh;
    o.mops.atrel = true;
    o.mops.atnorm = at_norminf;
    o.mops.atnodiag = false;


    FMatrix<double> hf(r, c, h);
//     RMatrix hs(hf, *comm);
    RMatrix hs(hf); // defaults to MPI_COMM_WORLD

    bool scr = false, sfg = false;

    LSChebHelper lsh;
    lsh.setops(o);
    lsh.setham(hs);
    lsh.compute(scr, sfg);

    FMatrix<double> rf(lsh.RHO);
    for (int i = 0; i < r*c; ++i) rho[i] = rf.access(i);
}




