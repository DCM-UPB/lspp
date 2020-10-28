#include "libls.hpp"
#include <fstream>

void LSHelper::setops(const LSOptions& nops)
{
    std::pout<<" * SETTING OPTIONS\n";
    pops=nops; 
    if (HAM.rows()>0) init_all();
}

void LSHelper::setham(const RMatrix& nham)
{
    std::pout<<" * SETTING THE HAMILTONIAN\n";
    HAM=nham;
    init_all(); 
}

void LSHelper::init_all() 
{
    std::pout<<" * COMPUTING INITIAL MATRICES\n";

    if (pops.flags & lsf_autohamspec) 
    {
        std::pout<<" ** EXTIMATING HAMILTONIAN SPECTRUM BOUND\n";
        /*RMatrix tmp,tmp2;
        MatrixOptions<RMatrix> hops; HAM.getops(hops); tmp.setops(hops); tmp2.setops(hops);
         std::cerr<<"now we compute the product\n"; 
        
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
    std::pout <<" ** COMPUTING EXPONENTIALS\n";
    RMatrix tmp(HAM); tmp-=pops.mu;
    
    MatrixOptions<CMatrix> mops;
    HAM.getops(mops);
    ETHAMP.setops(mops);
    tmp*=(-pops.beta/(2*pops.p));
    exp(tmp,ETHAMP, pops.accu*pops.accu);   //computes exponential with improved accuracy: it is cheap, so better being safe than sorry 
    
    finit=true;
}
/*
void LSHelper::seriescoeff_r(const double &q, const complex& k, std::valarray<double>& sc, const double& scale)
{
    unsigned long m=sc.size()-1;
    double cq=cos(constant::pi*(q+q-1)/(pops.p+pops.p)), sq=sqrt(1-cq*cq);
    complex psiq=complex(cq,sq);
    complex pk=complex(1,0)/(1.-k*psiq);
    complex ykq=pk*psiq*scale;  //scaling factor to avoid overflows...
    
    double bin;
    complex yi=complex(1,0),yj,u=complex(0,-k.imag()/scale),uj;  //scaling factor to avoid overflows...
    sc=0.;
    for (long i=0;i<=m;++i)
    {
        bin=1;  yj=yi; uj=1.;
        for (long j=i;j<=m;++j)
        { sc[i]+=real(pk*yj*bin*uj);  uj*=u; bin*=(j+1.)/(j+1.-i); yj*=ykq; }
        yi*=ykq;
    }
}
void LSHelper::tailcoeff_r(const double &q, const complex& k, std::valarray<double>& tc,const double& scale)
{
    unsigned long m=tc.size()-1;
    tc=0.;
    std::valarray<double> sc(tc);
    for(unsigned long qq=(unsigned long) q; qq<=(unsigned long) pops.p; ++qq) { seriescoeff_r(qq,  k, sc, scale); tc+=sc; }
}
*/
void LSHelper::seriescoeff(const double &q, const stdcomplex& k, std::valarray<stdcomplex>& sc,const double& scale)
{
    unsigned long m=sc.size()-1;
    double cq=cos(constant::pi*(q+q-1)/(pops.p+pops.p)), sq=sqrt(1-cq*cq);
    stdcomplex psiq=stdcomplex(cq,sq);
    stdcomplex yqk=stdcomplex(1,0)/(stdcomplex(cq,-sq)-k)*scale;  //scaling factor to avoid overflows
    stdcomplex pre=stdcomplex(1,0)/(1.-k*psiq);
    
    sc[0]=1.; if(m==0) return;
    sc[1]=yqk; for (unsigned long i=2;i<=m; ++i) {sc[i]=sc[i-1]; sc[i]*=yqk;}
    sc*=pre;
}
void LSHelper::tailcoeff(const double &q, const stdcomplex& k, std::valarray<stdcomplex>& tc,const double& scale)
{
    tc=stdcomplex(0.,0.);
    std::valarray<stdcomplex> sc(tc);
    for(unsigned long qq=(unsigned long) q; qq<=(unsigned long) pops.p; ++qq) { seriescoeff(qq,  k, sc, scale); tc+=sc; }
}


double re(const stdcomplex& c)
{
    return std::real(c);
}

void LSHelper::compute(bool& scratch, bool& safeguess) //mod by Doro
{
    if (!finit) init_all();
#ifdef BENCHMARK
    TBBenchmarks["LL_OVERALL"].n_calls++; TBBenchmarks["LL_OVERALL"].timer.start();
#endif
    
    // ******************** does the series part **********************  //
    //gets best k
    stdcomplex k; double ksv; 
    opt_k(pops.qbar, k, ksv);
    unsigned long ms=(unsigned long) ceil(log(pops.accu)/log(ksv));
    //plug here to do any test to check whether it is better to use std polynomial method
    bool ffastpoly=(pops.flags&lsf_fastpoly);
    if (ffastpoly && 3*sqrt(1. * ms) >= ms) 
    {
        ffastpoly=false;
        //std::cerr<<"Warning: fallback on standard polynomial evaluation, which should be faster here!\n";
    }
    std::cerr<<"Optimized k-1: "<<(k-1.)<<"\n";
    
    //we define both, then depending on options, different variables will be used.
    //I don't think there is significant overhead in choosing the mode at runtime
    //tail coefficients
    std::valarray<double> rtc(ms+1), rsc(ms+1);
    std::valarray<stdcomplex> ctc(ms+1), csc(ms+1);
    
    CMatrix cIMQ, cMQ, cRHO, cYk;
    cRHO.resize(HAM.rows(),HAM.rows()); 
    cIMQ.resize(HAM.rows(),HAM.rows()); 
    MatrixOptions<CMatrix> rmops;
    HAM.getops(rmops); 
    cRHO.setops(rmops); cIMQ.setops(rmops);
    cYk.setops(rmops); cMQ.setops(rmops);

    
    //computes Yk
    //to avoid overflow, we Yk is scaled by 1/|Yk|,
    //and the coefficients are scaled accordingly
    cYk=ETHAMP; cYk-=k;
    double nYk=norminf(cYk); cYk*=1./nYk;
    tailcoeff(pops.qbar, k, ctc, nYk);
    seriescoeff(pops.qbar,k,csc, nYk);
    
    //!RHO.resize(HAM.rows(), HAM.cols()); cIMQ.resize(0,0);  goto nwt;  //! temporary to check cheby
    std::cerr<<" * COMPUTING SERIES CONTRIBUTION, POLYNOMIAL ORDER: "<<ms<<"\n";
    if (pops.flags & lsf_niwarmstart)
    {
        std::valarray<std::valarray<stdcomplex>* > valc(2); 
        std::valarray<CMatrix *> valm(2); 
        valc[0]=&ctc; valc[1]=&csc;
        valm[0]=&cRHO;  
        valm[1]=&cIMQ;
        
        if (ffastpoly) poly_liang(cYk,valc,valm);
        else poly_standard_nonhorner(cYk,valc,valm);
    }
    else
    {
        if (ffastpoly) poly_liang(cYk,ctc,cRHO);
        else poly_standard_nonhorner(cYk,ctc,cRHO);
        cIMQ.resize(0,0); 
    }
    
/*    std::ofstream os;
    os.open("cRHO"); 
    os<<cRHO;
    os.close();
    
    os.open("rRHO"); 
    os<<rRHO;
    os.close();*/
    
    //we start storing the tail contrib.
    copymap(cRHO,RHO,re); cRHO*=0.;
    std::cerr<<" * SERIES CONTRIBUTION IS OVER, PARTIAL TRACE: "<<trace(RHO)/pops.p<<"\n";
    //to the newton part
    /* 
            newton inversion bailout options
      we consider inversion converged checking the norm of the residual:
      this NEEDs to be below 1
      if it is below accuracy OR if the change in the norm is below it, we quit
    */

    IterOptions<double,3> inwops(100,pops.accu, 0., ichk_default);
    inwops.thresh[0]=1.; //MINIMUM we want the residual to be below 1!
    inwops.flags[1]=ichk_sufficient;
    inwops.flags[2]=ichk_change;
    
    std::cerr<<" * COMPUTING NEWTON CONTRIBUTION\n";

    double cq=cos(constant::pi/pops.p), sq=sqrt(1-cq*cq);
    stdcomplex kextra=stdcomplex(cq,sq);
    CMatrix cIMQh1; cIMQh1.setops(rmops); cIMQh1.resize(0,0);
    
    for (unsigned long q=pops.qbar-1; q>0; --q)
    { 
        std::cerr<<" >> computing q-channel: "<<q<<"\n";
        cq=cos(constant::pi*(q+q-1)/(pops.p+pops.p)); sq=sqrt(1-cq*cq);
        stdcomplex psiq=stdcomplex(cq,sq);
        
        cMQ=ETHAMP; cMQ*=(-psiq); cMQ+=1;
        if (pops.flags & lsf_niextrapolation) {
            if (pops.flags & lsf_niextra2)
            {
                //long history extrapolation
                if (cIMQh1.rows()==0) cIMQ*=kextra; 
                else 
                {
                    cIMQ*=(kextra+kextra*kextra);
                    cIMQh1*=(-kextra*kextra*kextra);
                    cIMQ+=cIMQh1;
                }
                cIMQh1=cIMQ;
            }
            else
            {
                std::cerr<<"extrapolating, order 1\n";
                cIMQ*=kextra; 
            }
        } 
        else 
        {
            cIMQ.resize(0,0);
        }
        inverse_newton(cMQ,cIMQ,inwops);
        
        /*!
        std::string fname;
        std::stringstream ss;
        ss<<"invmat-"<<q;
        ss>>fname;
        std::ofstream os(fname.c_str());
        os<< cIMQ;
        os.close();
        */
        cRHO+=cIMQ;
        std::cerr<<"q: "<<q<<" tr(m_q^-1): "<<trace(cIMQ)/(1.*pops.p)<<"\n";
    }
     
    //stitches the pieces together
    RMatrix rRHO; copymap(cRHO,rRHO,re);
    RHO+=rRHO;
    RHO*=-(1./pops.p);
    RHO+=1.;
#ifdef BENCHMARK
    TBBenchmarks["LL_OVERALL"].timer.stop(); TBBenchmarks["LL_OVERALL"].tot_time+=TBBenchmarks["LL_OVERALL"].timer;
#endif

}

void LSHelper::opt_k(const double& q, stdcomplex& k, double &csi)
{
    fscratio fc;
    fc.set_pars(pops.p,q,(pops.hs_max-pops.mu)*pops.beta,(pops.hs_min-pops.mu)*pops.beta); 
    std::cerr<<"fcpars: q: "<< q<<" max: "<<fc.hs_max<<" min "<<fc.hs_min<<"\n";
    std::valarray<double> rp(2); double rv;
    //extimates based on theoretical arguments
    double hsr=fabs(pops.hs_min-pops.mu);
    if (hsr<fabs(pops.hs_max-pops.mu)) hsr=fabs(pops.hs_max-pops.mu); hsr*=pops.beta;
    rp[0]=(hsr*hsr)/(2*constant::pi*(q+q-1)); 
    rp[1]=constant::pi/2.;
    fc.set_vars(rp); fc.get_value(rv);
    
    if (pops.flags & lsf_optk)
    {
        std::cerr<<"optimizin' k; initial conv val: "<<rv<<" for hsr: "<<hsr<<"\n";
        min_simplex(fc,make_simplex(rp,0.1,1.1),rp,rv, 
                IterOptions<double,2>(
                                      10000,
                                      fixarray<double,2>(__TB_STD_EPS, __TB_STD_EPS),
                                      fixarray<double,2>(0.,0.),
                                      fixarray<double,2>(ichk_change, ichk_default))
                   );
        std::cerr<<"optimizin' k; final conv val: "<<rv<<"\n";
    }
    
    csi=sqrt(rv); k=stdcomplex(1+rp[0]/pops.p*cos(rp[1]), rp[0]/pops.p*sin(rp[1]));
}

void LSHelper::opt_qbar()
{
    if (pops.flags & lsf_niwarmstart) pops.qbar=2; else pops.qbar=1;
    double nm_series, nm_newt;
    double hsr=fabs(pops.hs_min-pops.mu);
    if (hsr<fabs(pops.hs_max-pops.mu)) hsr=fabs(pops.hs_max-pops.mu); hsr*=pops.beta;
    while (pops.qbar<pops.p)
    {
        //compute conv. speed for series
        stdcomplex ok; double ov, ov_next;
        opt_k(pops.qbar, ok, ov);
        opt_k(pops.qbar+1, ok, ov_next);
        
        //compute conv. for newton
        double cn, nm, hm;
        if (pops.flags & lsf_niwarmstart)
        {
            hm=2./(2*pops.qbar-1);
        }
        else
        {
            
            //cold start extimates
            cn=1+hsr/(constant::pi*(pops.qbar+pops.qbar-1)); //condition number of m(q)
            hm=1.-1./(cn*cn*HAM.rows());
        }
        //extimate of the spectral radius of m_q
        nm=sqrt(hsr*hsr+constant::pi*constant::pi*(pops.qbar+pops.qbar-1)*(pops.qbar+pops.qbar-1))/(pops.p+pops.p);
        nm_newt=2*ceil(1/log(2.)*log(log(pops.accu*nm*(1-hm))/log(hm)));
        
        //if you consider the fast poly...
        if (pops.flags & lsf_fastpoly)
        {
            nm_series=ceil(sqrt(log(pops.accu)/log(ov)))-ceil(sqrt(log(pops.accu)/log(ov_next)));
            if (pops.flags & lsf_niwarmstart) nm_series*=3.; else nm_series*=2.;
        }
        else
            nm_series=ceil(log(pops.accu)*(1/log(ov)-1/log(ov_next)));
        //if (pops.flags & lsf_fastpoly) nm_series=ceil(sqrt(nm_series))*(pops.flags & lsf_niwarmstart?3:2);
        
        std::cerr<<"Polynomial order: " << log(pops.accu)*(1/log(ov))<<"\n";
        std::cerr<<"For qbar= "<<pops.qbar<<" estimated n. newton: "<<nm_newt<<",  n. series: "<<nm_series<<"\n";
        if (nm_series<nm_newt) break;
        pops.qbar++;
    }
}

