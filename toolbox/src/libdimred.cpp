#include "dimreduce.hpp"
#include "matrix-io.hpp"
#include "matrix-conv.hpp"
#include "linalg.hpp"

#ifdef ARPACK
namespace tblapack {
extern "C" {
#define dsaupd dsaupd_
#define dseupd dseupd_
    int dsaupd( int* IDO, char* BMAT, int*N, char* WHICH, int* NEV, 
             double* TOL, double *RESID, int* NCV, double* V, int *LDV, 
             int* IPARAM, int* IPNTR, double *WORKD, 
             double *WORKL, int*LWORKL, int*INFO );
    int dseupd( int* RVEC, char* HOWMNY, int *SELECT, 
                double *D, double *Z, int *LDZ, double *SIGMA, 
                char* BMAT, int*N, char* WHICH, int* NEV, 
                double* TOL, double *RESID, int* NCV, double* V, int *LDV, 
                int* IPARAM, int* IPNTR, double *WORKD, 
                double *WORKL, int*LWORKL, int*INFO);
    
}};
#endif

namespace toolbox {
/*! A collection of metric functions */
double NLDRMetricEuclid::pdist(const double* a, const double* b, unsigned long n) const
{
    double d=0.0;
    for (unsigned long i=0; i<n; ++i) d+=(b[i]-a[i])*(b[i]-a[i]);
    return sqrt(d);
}

void NLDRMetricPBC::pdiff(const double* a, const double* b, double* c, unsigned long n) const
{
#ifdef DEBUG
    if (n!=periods.size()) ERROR("Periodicity array has wrong dimensions\n");
#endif
    double dx;

    for (unsigned long i=0; i<n; ++i) 
    { 
        dx=(b[i]-a[i]); 
        dx/=periods[i]; 
        dx-=round(dx);  
        dx*=periods[i]; 
        c[i]=dx;
    }
}
double NLDRMetricPBC::pdist(const double* a, const double* b, unsigned long n) const
{
#ifdef DEBUG
    if (n!=periods.size()) ERROR("Periodicity array has wrong dimensions\n");
#endif
    double d=0.0, dx;

    for (unsigned long i=0; i<n; ++i) 
    { 
        dx=(b[i]-a[i]);
        dx/=periods[i];
        dx-=round(dx);
        dx*=periods[i];
        d+=dx*dx;
    }
    return sqrt(d);
}

/*Builds neighbor list*/
void NLDRNeighborList::nlbuildup(const FMatrix<double>& points)
{
    unsigned long n=points.rows(), d=points.cols(), tn=0, pn;
    std::valarray<NLDRNeighbor> ni(n); 
    std::valarray<std::vector<NLDRNeighbor> > nn(n);
    
    npoint.resize(n+1); 
    
    //we rather compute distances twice but avoid storing a matrix which is n x n
    for (unsigned long i=0; i<n; i++)
    {
        //builds array with all the distances
        for (unsigned long j=0; j<n; j++) {
            ni[j].j=j; 
            ni[j].d=opts.ometric->pdist(&(const_cast<FMatrix<double>&>(points)(i,0)), &(const_cast<FMatrix<double>&>(points)(j,0)), d); }
        heapsort(ni); //sorts distances
        //picks maxneigh and/or within cutoff neighbors (skips self)
        unsigned long k=0; 
        while (
               ((opts.kw>0 && k<opts.kw) || ( 
               (opts.cutoff==0.0 || ni[k+1].d<opts.cutoff ) &&  
                 (opts.maxneigh==0 || (opts.kw==0 && k<opts.maxneigh)) ))
                && k<n) ++k;
        
        nn[i].resize(k); for (unsigned long j=0; j<k; ++j) nn[i][j]=ni[j+1];
        npoint[i]=k;  tn+=k; 
    }
    if (opts.kw>0)
    {
        std::cerr<<"Rebuilding list with weighted distances.\n";
        //variables for WLLE metric
        double wllec1, wllec2; double wllea, wlleb, wllel; std::valarray<double> wlletau(d);
        wllec2=sqrt(2.0)*exp(lgamma((d+1)*0.5)-lgamma(d*0.5));
        wllec1=wllec2/d; tn=0;
        std::cerr<<"WLLEC1,2  "<<wllec1<<","<<wllec2<<"\n";
        //builds deformed distance à la WLLE - which needs cycling again over the neighbors
        for (unsigned long i=0; i<n; i++)
        {
            //computes for the metric
            wlletau=0.0; wllel=0.0;
            for (unsigned long j=0; j<nn[i].size(); ++j)
            { 
                for (unsigned long k=0; k<d; ++k) wlletau[k]+=points(nn[i][j].j,k); 
                wllel+=nn[i][j].d;
            }
            wlletau*=1.0/nn[i].size(); for (unsigned long k=0; k<d; ++k) wlletau[k]-=points(i,k);
            wllel*=1.0/nn[i].size();
            
            wllea=wllel/wllec2; 
            wlleb=0.0; for (unsigned long k=0; k<d; ++k) wlleb+=wlletau[k]*wlletau[k]; wlleb=sqrt(wlleb);
            //std::cerr<<"L,  |avg| "<< wllel<<" "<<wlleb<<"\n";
            wlletau*=1.0/wlleb; wlleb*=1.0/wllec1;
            //std::cerr<<"WLLE a,b "<<wllea<<","<<wlleb<<"\n";
            if (wlleb>wllea) wlleb=wllea; //avoids getting negative distances just because the data are not CAM-distributed! this method just sucks.
            //builds array with all the distances
            for (unsigned long j=0; j<n; j++) {
                
                ni[j].j=j; 
                if (i==j) { ni[j]=0.0; continue; }
                ni[j].d=opts.ometric->pdist(&(const_cast<FMatrix<double>&>(points)(i,0)), &(const_cast<FMatrix<double>&>(points)(j,0)), d); 
                //std::cerr<<i<<","<<j<<" "<<ni[j].d<<" >> ";
                double wlletx=0.0; for (unsigned long k=0; k<d; ++k) wlletx+=(points(j,k)-points(i,k))*wlletau[k];
                ni[j].d=ni[j].d/(wllea+wlleb*wlletx/ni[j].d);
                //std::cerr<<ni[j].d<<"\n";
            }

            heapsort(ni); //sorts distances
            //picks maxneigh and/or within cutoff neighbors (skips self)
            unsigned long k=0; 
            while ( (opts.cutoff==0.0 || ni[k+1].d<opts.cutoff ) &&  
                     (opts.maxneigh==0 || k<opts.maxneigh ) && k<n) ++k;
    
            nn[i].resize(k); for (unsigned long j=0; j<k; ++j) nn[i][j]=ni[j+1];
            npoint[i]=k;  tn+=k; 
        }
    }
        
    //implements greedy/liberal symmetrization of the neighbor list
    if (opts.greediness!=NLDRAsym)
    for (unsigned long i=0; i<n; ++i)
    {
        for (unsigned long j=0; j<npoint[i]; ++j)
        {
            bool ijsym=false; int ijj=nn[i][j].j;
            for (unsigned long k=0; k<npoint[ijj]; ++k) if (nn[ijj][k].j==i) { ijsym=true; break; } 
            if (!ijsym) {
                switch(opts.greediness) {
                    case NLDRGreedy:
                        nn[ijj].push_back(NLDRNeighbor(i,nn[i][j].d)); 
                        ++npoint[ijj]; ++tn;
                        break;
                    case NLDRLiberal:
                        nn[i].erase(nn[i].begin()+j);
                        --npoint[i]; --j; --tn;
                        break;
                    default:
                        ERROR("Unsupported NL symmetrization scheme");
                }
            }
        }
    }
    
    //collapses the (expensive) storage as vectors to a "compressed" storage
    nlist.resize(tn); unsigned long k=0; pn=tn=0;
    for (unsigned long i=0; i<n; ++i)
    {
        pn=tn; tn+=npoint[i]; std::valarray<NLDRNeighbor> nni(nn[i].size()); 
        for (unsigned long j=0; j<npoint[i]; ++j) nni[j]=nn[i][j]; if (npoint[i]>1) heapsort(nni);
        for (unsigned long j=0; j<npoint[i]; ++j) nlist[k++]=nni[j];
        npoint[i]=pn;
    }
    npoint[n]=tn;
}

NLDRNeighbor& NLDRNeighborList::rneigh(unsigned long i, unsigned long j)
{
#ifdef DEBUG
    if (i>npoint.size()-1 || j>(npoint[i+1]-npoint[i])) ERROR("Trying to access non-existent neighbor.");
#endif 
    return nlist[npoint[i]+j];
}

unsigned long NLDRNeighborList::index(unsigned long i, unsigned long j) const 
{ return const_cast<NLDRNeighborList*>(this)->rneigh(i,j).j; }
double NLDRNeighborList::dist(unsigned long i, unsigned long j) const
{ return const_cast<NLDRNeighborList*>(this)->rneigh(i,j).d; }
NLDRNeighbor NLDRNeighborList::neigh(unsigned long i, unsigned long j) const
{ return const_cast<NLDRNeighborList*>(this)->rneigh(i,j); }
void NLDRNeighborList::RemoveLonesome(std::valarray<unsigned long>& llone)
{
    unsigned long nlone=0;
    for (unsigned long i=0; i<npoint.size()-1; ++i) if ((npoint[i+1]-npoint[i])<=0) ++nlone;
    llone.resize(nlone); if (nlone==0) return;
    std::valarray<unsigned long> newnp(npoint.size()-nlone);
    nlone=0;
    for (unsigned long i=0; i<npoint.size()-1; ++i) 
    {
        newnp[i-nlone]=npoint[i];
        if ((npoint[i+1]-npoint[i])<=0) { llone[nlone]=i; ++nlone; }
    }
    
    //then updates indices
    for (unsigned long i=0; i<nlist.size(); ++i) 
    { 
        unsigned long dni=0;
        for (unsigned long j=0; j<nlone; ++j) if (nlist[i].j>llone[j]) ++dni;  else break;
        nlist[i].j-=dni;
    }
    
    newnp[newnp.size()-1]=npoint[npoint.size()-1];
    npoint.resize(newnp.size()); npoint=newnp;
}

double isqrt(double x) { return 1.0/sqrt(x); }
void NLDRLLE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLEOptions& opts, NLDRLLEReport& report)
{
    unsigned long dts=(opts.dimts==0?opts.lowdim:opts.dimts);
    proj.P=points; proj.D=points.cols(); proj.d=opts.lowdim; proj.n=points.rows();
    proj.p.resize(proj.n,proj.d);
    
    std::cerr<<"Building neighbor list\n";
    NLDRNeighborList nlist(proj.P,opts.nlopts);
    if (opts.rmlonesome) 
    {
        std::valarray<unsigned long> ilone;
        nlist.RemoveLonesome(ilone);
        //removing points with no connections
        unsigned long nlone=ilone.size();
            
        if (nlone>0) 
        {
            std::cerr<<"Removing isolated points\n";
            proj.P.resize(proj.n-nlone,proj.D);
            unsigned long k=0, j=0;
            for (unsigned long i=0; i<proj.n; ++i) 
            {
                if (j<nlone && ilone[j]==i) { ++j; continue; }
                for (unsigned long h=0; h<proj.D; ++h) proj.P(k,h)=points(i,h);
                ++k;
            }
            proj.n-=nlone;
            proj.p.resize(proj.n,proj.d);
        }
    }
    
    
    std::cerr<<"Building weights\n";
    std::valarray<std::valarray<double> > weights(proj.n);
    std::valarray<std::valarray<unsigned long> > indices(proj.n);
    
    if (opts.verbose) { report.hd_errors.resize(proj.n); report.hd_errors=0.0;  report.ld_errors.resize(proj.n); report.ld_errors=0.0; }
    report.ld_error=report.hd_error=0.0;
    
    for (unsigned long i=0; i<proj.n; ++i)  {  unsigned long m=nlist.nneigh(i); if (m==0) ERROR("Point "<<i<<" is isolated!\n"); weights[i].resize(m); indices[i].resize(m); }
    
    CrsMatrix<double> W, WT; double dp=dts*(dts+1)/2; std::valarray<double> x1(proj.D), x2(proj.D), dx12(proj.D);
    if (opts.mode==HLLE) W.resize(proj.n*dp,proj.n);
    for (unsigned long i=0; i<proj.n; ++i)
    {
        unsigned long m=nlist.nneigh(i);

        FMatrix<double> C(m,m), C1, G(proj.D,m), GT, H, HT;
        x1.resize(proj.D); x2.resize(proj.D); 
        x1=proj.P.row(i);
        //builds matrix with neighbor vector shifts around the central point
        for (unsigned long j=0; j<m; ++j) 
        {
            x2=proj.P.row(nlist.index(i,j));
            opts.nlopts.ometric->diff(x2,x1,dx12);

            G.col(j)=dx12;
        } 
        transpose(G,GT);

        std::cerr<<"POINT "<<i<<": ";
        for (unsigned long j=0; j<m; ++j)  std::cerr<<nlist.index(i,j)<<"("<<nlist.dist(i,j)<<") ";
        std::cerr<<"\n";
            
        if(opts.mode==LLE)
        {
             mult(GT,G,C);
        }
        else if (opts.mode==LLTE)
        {
            //First, finds d principal components
            mult(G,GT,C);
            FMatrix<double> P, PT; std::valarray<double> p;
            EigenSolverSym(C,P,p,LAEMIndex,proj.D-dts,proj.D-1);

            transpose(P,PT); mult(PT,G,H);  transpose(H, HT);
            mult(HT,H,C);  // now C is projected on the most relevant eigenvalues
        }
        else if (opts.mode==HLLE)
        {
            mult(GT,G,C);
            FMatrix<double> U; std::valarray<double> u;
            EigenSolverSym(C,U,u,LAEMIndex,m-dts,m-1);
            
            //builds Hessian estimator
            x1.resize(m); x2.resize(m);
            FMatrix<double> Xi(m,1+dts+dp);
            x1=1.; Xi.col(0)=x1;
            for (unsigned long h=0; h<dts; ++h) Xi.col(1+h)=U.col(h);
            unsigned long k=0;
            
            for (unsigned long h=0; h<dts; ++h) { x1=U.col(h); for (unsigned long h2=h; h2<dts; ++h2) 
            {
                x2=U.col(h2); x2*=x1;
                Xi.col(1+dts+k)=x2; k++;
            } }

            //GS-orthogonalize columns
            FMatrix<double> R(Xi.cols(),Xi.cols(),0.0);
            for (unsigned long h1=0; h1<R.rows(); h1++)
            {
                R(h1,h1)=0.0; for (unsigned long h2=0; h2<Xi.rows(); h2++) R(h1,h1)+=Xi(h2,h1)*Xi(h2,h1); R(h1,h1)=sqrt(R(h1,h1));
                for (unsigned long h2=0; h2<Xi.rows(); h2++) Xi(h2,h1)*=1./R(h1,h1);
                if (h1==R.rows()-1) break;
                for (unsigned long h2=h1+1; h2<R.rows(); h2++) 
                {
                    R(h1,h2)=0.0; for (unsigned long h3=0; h3<Xi.rows(); h3++) R(h1,h2)+=Xi(h3,h1)*Xi(h3,h2);
                    for (unsigned long h3=0; h3<Xi.rows(); h3++) Xi(h3,h2)-=R(h1,h2)*Xi(h3,h1);
                }
            }

            H.resize(dp,m);
            for (unsigned long h=0; h<dp; ++h) H.row(h)=Xi.col(h+1+dts);
        }
        else ERROR("Unsupported LLE mode");

        if (opts.mode==LLE || opts.mode==LLTE) 
        {
            //smoothens the matrix
            double tr=0.0; 
            if (opts.smooth<0) { tr=normfrob(C)*(-opts.smooth); }
            else tr=opts.smooth;
            for (unsigned long j=0; j<m; ++j) C(j,j)+=tr;
            
            MatrixInverse(C,C1);
            
            double beta, lambda;
            beta=0.0;
            for (unsigned long j=0; j<m; ++j) for (unsigned long k=0; k<m; ++k) beta+=C1(j,k);
            lambda=1.0/beta;
            for (unsigned long j=0; j<m; ++j)
            {
                indices[i][j]=nlist.index(i,j);  weights[i][j]=0.0; 
                for (unsigned long k=0; k<m; ++k) weights[i][j]+=C1(j,k);
            }
            weights[i]*=lambda;
            
            //with LLTE we are computing the error in the reconstruction of PROJECTED high-D points
            for (unsigned long h=0; h<(opts.mode==LLE?proj.D:proj.d); ++h) 
            { 
                double ed=0.0;
                if (opts.mode==LLE)
                {
                    for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*proj.P(nlist.index(i,j),h);
                    ed-=proj.P(i,h); ed*=ed;
                }
                else if (opts.mode==LLTE)
                {
                    for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*H(h,j);
                    ed*=ed;  
                }
                report.hd_error+=ed; if (opts.verbose) report.hd_errors[i]+=ed;
            }
            if (opts.verbose) report.hd_errors[i]=sqrt(report.hd_errors[i]);
        }
        else if (opts.mode==HLLE)
        {
            if (opts.verbose)
            {
                FMatrix<double> HT, HHT;
                transpose(H,HT); mult(HT,H,HHT);
                report.hd_errors[i]=trace(HHT)/HHT.rows();
                FMatrix<double> U; std::valarray<double> u;
                EigenSolverSym(HHT,U,u);
                std::cerr<<" singular values of HHT for point "<<i<<"\n"<<u;
                std::cerr<<report.hd_errors[i]<<" "<<normfrob(H)<<"\n";
            }
/*            FMatrix<double> HT, HHT, U; std::valarray<double> u;
            transpose(H,HT); mult(HT,H,HHT);
            EigenSolverSym(C,U,u);
            report.hd_errors[i]=u[u.size()-1];
            std::cerr<<" singular values of HHT for point "<<i<<"\n";
            //if (i==1130 || i==833) H*=0.0; 
            std::cerr<<u;
*/
            //builds W straight away
            for (unsigned long h=0; h<dp; ++h)
            {
                for (unsigned long j=0; j<m; ++j)
                    W(i*dp+h,nlist.index(i,j))=H(h,j);
            }
        }
    }
    report.hd_error=sqrt(report.hd_error/proj.n);
    
    
    if (opts.mode==LLE||opts.mode==LLTE)
    {
        //we have the weight matrix as a sparse matrix
        std::cerr<<"Building sparse W matrix\n";
        W=CrsMatrix<double>(proj.n,proj.n,indices,weights);
        W*=-1.0;  for (unsigned long i=0; i<proj.n; ++i) W(i,i)+=1.0;
    }
    //finds eigenvalues by calling ARPACK routines, else uses full matrix algebra...
#ifdef ARPACK
    int IDO, N, NEV, NCV, LDV, IPARAM[11], IPNTR[11], LWORKL, INFO;
    char BMAT, WHICH[2]; 
    BMAT='I'; WHICH[0]='S'; WHICH[1]='A';
    NEV=proj.d+2; N=proj.n;
    NCV=50*(NEV+1); if(NCV>N) NCV=N; 
    LDV=N;
    for (int i=0; i<11; i++) IPARAM[i]=IPNTR[i]=0;
    IPARAM[0]=1; IPARAM[2]=3*N; IPARAM[3]=1; 
    IPARAM[4]=NEV; IPARAM[6]=1; 
    double TOL=1e-5; 
    IDO=0; LWORKL=NCV*(NCV+8); INFO=0; 
    std::valarray<double> V(N*NCV), RESID(N), WORKD(3*N), WORKL(LWORKL); 
    V=0.0; RESID=0.0; WORKD=0.0; WORKL=0.0;
    std::valarray<double> y(N);
    IDO=0;
    
    do {
/*        std::cerr<<"Calling DSAUPD\n";
        std::cerr<<"IDO "<<IDO<<"\n";
        std::cerr<<"N "<<N<<"\n";
        std::cerr<<"NEV "<<NEV<<"\n";
        std::cerr<<"NCV "<<NCV<<"\n";
        std::cerr<<"LDV "<<LDV<<"\n";
        std::cerr<<"LWORKL "<<LWORKL<<"\n";
        std::cerr<<"INFO "<<INFO<<"\n";
        std::cerr<<"IPARAM "; for (int i=0; i<11; i++)std::cerr<<IPARAM[i]<<" "; std::cerr<<"\n";
        std::cerr<<"IPNTR "; for (int i=0; i<11; i++)std::cerr<<IPNTR[i]<<" "; std::cerr<<"\n";
        */
        tblapack::dsaupd(&IDO, &BMAT, &N, &WHICH[0], &NEV, 
               &TOL, &RESID[0], &NCV, &V[0], &LDV, 
               &IPARAM[0], &IPNTR[0], &WORKD[0], 
               &WORKL[0], &LWORKL, &INFO);
        
        std::valarray<double> x(WORKD[std::slice(IPNTR[0]-1,N,1)]); //takes slice
        mult(W,x,y); x=y; Tmult(W,x,y);
        std::cerr<<"Returned IDO: "<<IDO<<" INFO: "<<INFO<<"\n";
        if (INFO!=0) std::cerr<<" IPARAM[4]: "<<IPARAM[4]<<"\n";
        WORKD[std::slice(IPNTR[1]-1,N,1)]=y;
    } while (IDO==1 || IDO==-1);
    
    int RVEC; std::valarray<int> SELECT(NCV); std::valarray<double> D(N), Z(N*NEV);
    int LDZ; char HOWMNY; double SIGMA;
    RVEC=1; HOWMNY='A'; LDZ=N; SIGMA=0.0; SELECT=0; D=0.0; Z=0.0; 
    
    tblapack::dseupd(&RVEC, &HOWMNY, &SELECT[0], &D[0], &Z[0], &LDZ, &SIGMA,
            &BMAT, &N, WHICH, &NEV, 
            &TOL, &RESID[0], &NCV, &V[0], &LDV, 
            IPARAM, IPNTR, &WORKD[0], 
            &WORKL[0], &LWORKL, &INFO);
    for (int i=0; i<NEV; i++) std::cerr<<"EIGV: "<<i<<" is "<<D[i]<<"\n";
    
    for (unsigned long i=0; i<proj.n; ++i) for (unsigned long h=0; h<proj.d; ++h) 
        proj.p(i,h)=Z[i+(h+1)*proj.n];
#else
    std::cerr<<"Finding low-dim points (matrix-matrix mult)\n";
    transpose(W,WT); CrsMatrix<double> MM; mult(WT,W,MM); 
    FMatrix<double> FM(MM), Q; std::valarray<double> q;
    
    report.ld_errors=FM.diag();
    
    std::valarray<double> c(proj.n);
    for (unsigned long i=0; i<proj.n; ++i) if (FM(i,i)<0.5)
    {
        c=FM.col(i); c*=c; std::cerr<<"element "<<i<<" diag "<<FM(i,i)<<" row |sum| "<<sqrt(c.sum())<<"\n";
    }
        
    /*
    double dw=2;
    for (unsigned long i=0; i<proj.n; ++i) if (FM(i,i)<2)
    for (unsigned long j=0; j<proj.n; ++j) 
    { FM(i,j)-=dw/proj.n; if(i==j) FM(i,i)+=dw; else FM(j,i)-=dw/proj.n; }
    */
    
    /*FM.diag()+=wbad;
    wbad*=(1.0/proj.n);
    for (unsigned long i=0; i<proj.n; ++i) FM.row(i)-=wbad;*/
    //std::cerr<<"DIAGONAL: "<<wbad<<"\n";
    //FM(833,833)+=1e-3; FM(1130,1130)+=1e-3;
    /*FM.row(833)*=wbad; FM.col(833)*=wbad; 
    FM.row(1130)*=wbad; FM.col(1130)*=wbad; 
    std::cerr<<"Scaled elements of FM "<<FM(833,0)<<"; "<<FM(1130,1130)<<"\n";
    FM.row(100)*=wbad; FM.col(100)*=wbad; 
    FM.row(101)*=wbad; FM.col(101)*=wbad; 
    FM.row(102)*=wbad; FM.col(102)*=wbad; */
    
    //std::cerr<<FM<<"\n";
    std::cerr<<"Finding low-dim points (diagonalization "<<FM.rows()<<"x"<<FM.cols()<<")\n";
    if (opts.verbose) {
        EigenSolverSym(FM,Q,q,LAEMIndex,proj.d+1,proj.d+1);
        report.dp1eval=q[0];
    }
    //std::cerr<<"Finding zero eigenvec\n";
    //EigenSolverSym(FM,Q,q,LAEMIndex,0,0);
    //std::cerr<<"zero-eigenval: "<<q<<"eigenvec"<<Q<<"\n";
    EigenSolverSym(FM,Q,q,LAEMIndex,1,proj.d);
    //std::cerr<<"Finding zero eigenvec\n";
    report.deval.resize(proj.d); report.deval=q;
#endif
    
    if (opts.mode==LLE||opts.mode==LLTE)
    {
        for (unsigned long i=0; i<proj.n; ++i) for (unsigned long h=0; h<proj.d; ++h) proj.p(i,h)=Q(i,h);

        for (unsigned long i=0; i<proj.n; ++i)
        {
            unsigned long m=nlist.nneigh(i);
            for (unsigned long h=0; h<proj.d; ++h) 
            { 
                double ed=0.0;
                for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*proj.p(nlist.index(i,j),h);
                ed-=proj.p(i,h); ed*=ed;
                report.ld_error+=ed; if (opts.verbose) report.ld_errors[i]+=ed*ed;
            }
            if (opts.verbose) report.ld_errors[i]=sqrt(report.ld_errors[i]);
        }
    }
    else if (opts.mode==HLLE) 
    {
        std::cerr<<"HLLE Orienting coordinates\n";
        Q*=1.0/sqrt(proj.n);
        FMatrix<double> QT, R, R2; transpose(Q,QT); mult(QT,Q,R); MatrixFunctionSym(R,&isqrt,R2);
        mult(Q,R2,R);
        //R=Q; // no scaling
        for (unsigned long i=0; i<proj.n; ++i) for (unsigned long h=0; h<proj.d; ++h) proj.p(i,h)=R(i,h);
    }
}

// void NLDRLLTE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLTEOptions& opts, NLDRLLTEReport& report)
// {
//     proj.P=points; proj.D=points.cols(); proj.d=opts.lowdim; proj.n=points.rows();
//     proj.p.resize(proj.n,proj.d);
//     
//     std::cerr<<"Building neighbor list\n";
//     NLDRNeighborList nlist(proj.P,opts.nlopts);
//     //std::cerr<<proj.P<<"BEFORE \n";
//     if (opts.rmlonesome) 
//     {
//         std::valarray<unsigned long> ilone;
//         nlist.RemoveLonesome(ilone);
//         //removing points with no connections
//         unsigned long nlone=ilone.size();
//             
//         if (nlone>0) 
//         {
//             std::cerr<<"Removing isolated points\n";
//             proj.P.resize(proj.n-nlone,proj.D);
//             unsigned long k=0, j=0;
//             for (unsigned long i=0; i<proj.n; ++i) 
//             {
//                 if (j<nlone && ilone[j]==i) { ++j; continue; }
//                 for (unsigned long h=0; h<proj.D; ++h) proj.P(k,h)=points(i,h);
//                 ++k;
//             }
//             proj.n-=nlone;
//             proj.p.resize(proj.n,proj.d);
//         }
//     }
//     
//     std::cerr<<"Building weights\n";
//     std::valarray<std::valarray<double> > weights(proj.n);
//     std::valarray<std::valarray<unsigned long> > indices(proj.n);
//     
//     if (opts.verbose) { report.hd_errors.resize(proj.n); report.hd_errors=0.0;  report.ld_errors.resize(proj.n); report.ld_errors=0.0; }
//     report.ld_error=report.hd_error=0.0;
//     
//     for (unsigned long i=0; i<proj.n; ++i)
//     {
//         unsigned long m=nlist.nneigh(i);
//         if (m==0) ERROR("Point "<<i<<" is isolated!\n");
//         FMatrix<double> G(proj.D,m), GT, H, HT, C, C1; 
//         
//         for (unsigned long h=0; h<proj.D; ++h) for (unsigned long j=0; j<m; ++j) G(h,j)=points(nlist.index(i,j),h)-points(i,h);
//         
//         //First, finds d principal components
//         transpose(G,GT); mult(G,GT,C);
//         FMatrix<double> P, PT; std::valarray<double> p;
//         EigenSolverSym(C,P,p,LAEMIndex,proj.D-proj.d,proj.D-1);
// 
//         transpose(P,PT); mult(PT,G,H);  transpose(H, HT);
//         mult(HT,H,C);  // now C is projected on the most relevant eigenvalues
//         
//         //smoothens the matrix
//         double tr=0.0; 
//         if (opts.smooth<0) { tr=normfrob(C)*-opts.smooth; }
//         else tr=opts.smooth;
//         for (unsigned long j=0; j<m; ++j) C(j,j)+=tr;
//         
//         MatrixInverse(C,C1); 
//         
//         double beta, lambda;
//         beta=0.0;
//         for (unsigned long j=0; j<m; ++j) for (unsigned long k=0; k<m; ++k) beta+=C1(j,k);
//         lambda=1.0/beta;
// 
//         weights[i].resize(m); indices[i].resize(m);
//         for (unsigned long j=0; j<m; ++j)
//         {
//             indices[i][j]=nlist.index(i,j);  weights[i][j]=0.0; 
//             for (unsigned long k=0; k<m; ++k) weights[i][j]+=C1(j,k);
//         }
//         weights[i]*=lambda;
// 
//         for (unsigned long h=0; h<proj.d; ++h) 
//         { 
//             double ed=0.0;
//             for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*H(h,j);
//             ed*=ed;  report.hd_error+=ed; if (opts.verbose) report.hd_errors[i]+=ed*ed;
//         }
//         //should also compute TRUE HD errors, not just the projected ones!
//         if (opts.verbose) report.hd_errors[i]=sqrt(report.hd_errors[i]);
//     }
//     report.hd_error=sqrt(report.hd_error/proj.n);
//     
//     //we have the weight matrix as a sparse matrix
//     std::cerr<<"Building sparse W matrix\n";
//     CrsMatrix<double> W(proj.n,proj.n,indices,weights), WT;
//     W*=-1.0;  for (unsigned long i=0; i<proj.n; ++i) W(i,i)+=1.0;
// 
// 
//     std::cerr<<"Finding low-dim points\n";
//     transpose(W,WT); CrsMatrix<double> M; mult(WT,W,M); 
//     FMatrix<double> FM(M), Q; std::valarray<double> q;
//     EigenSolverSym(FM,Q,q,LAEMIndex,1,proj.d);
//     for (unsigned long i=0; i<proj.n; ++i) for (unsigned long h=0; h<proj.d; ++h) 
//             proj.p(i,h)=Q(i,h);
//     
//     report.deval.resize(proj.d); report.deval=q;
//     if (opts.verbose) {
//         EigenSolverSym(FM,Q,q,LAEMIndex,proj.d+1,proj.d+1);
//         report.dp1eval=q[0];
//     }
//     
//     for (unsigned long i=0; i<proj.n; ++i)
//     {
//         unsigned long m=nlist.nneigh(i);
//         for (unsigned long h=0; h<proj.d; ++h) 
//         { 
//             double ed=0.0;
//             for (unsigned long j=0; j<m; ++j) ed+=weights[i][j]*proj.p(nlist.index(i,j),h);
//             ed-=proj.p(i,h); ed*=ed;
//             report.ld_error+=ed; if (opts.verbose) report.ld_errors[i]+=ed*ed;
//         }
//         if (opts.verbose) report.ld_errors[i]=sqrt(report.ld_errors[i]);
//     }
// }

} //ends namespace toolbox