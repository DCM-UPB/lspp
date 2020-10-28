#include "dimreduce.hpp"
#include "clparser.hpp"
#include "matrix-io.hpp"

using namespace toolbox;
int main(int argc, char**argv)
{
    CLParser clp(argc,argv);
    unsigned long D,d,dts,nn; double sm, neps,peri; bool fverb, fveryverb;
    std::string nlrmode;
    bool fok=clp.getoption(D,"D",(unsigned long) 3) && 
            clp.getoption(d,"d",(unsigned long) 2) &&
            clp.getoption(dts,"dts",(unsigned long) 0) &&
            clp.getoption(sm,"smooth",-1e-3) &&  
            clp.getoption(nn,"neigh",(unsigned long) 4) &&
            clp.getoption(neps,"ncut",0.0) &&
            clp.getoption(peri,"period",0.0) &&
            clp.getoption(fverb,"v",false) &&  
            clp.getoption(nlrmode,"mode",std::string("LLE")) &&  
            clp.getoption(fveryverb,"vv",false);
    if (dts==0) dts=d; 
    std::vector<std::vector<double> > plist; std::vector<double> point(D);
    
    // reads points from standard input
    while (std::cin.good())
    {
        for (int i=0; i<D; i++) std::cin>>point[i];
        if (std::cin.good()) plist.push_back(point);
    }
    
    std::valarray<std::valarray<double> > hplist, lplist; 
    FMatrix<double> mpoints(plist.size(),D);
    for (int i=0; i<plist.size(); i++) for (int j=0; j<D; j++) mpoints(i,j)=plist[i][j];
    
    NLDRProjection nlproj;
    NLDRNeighborOptions nopts; nopts.greediness=NLDRAsym; nopts.maxneigh=nn; nopts.cutoff=neps; 
    NLDRMetricPBC nperi;
    if (peri==0.0) nopts.ometric=new NLDRMetricEuclid;
    else { 
        nperi.periods.resize(D); nperi.periods=peri;
        nopts.ometric=&nperi;  }
    NLDRLLEReport llereport;
    NLDRLLEOptions lleopts; lleopts.nlopts=nopts; lleopts.smooth=sm;  lleopts.lowdim=d; lleopts.dimts=dts;
    lleopts.rmlonesome=true; lleopts.verbose=fveryverb; 
    
    /*
    std::valarray<double> nr(3), nz(3); nr=1.0; nz=0.0;
    std::cerr<<lleopts.nlopts.ometric->dist(nr,nz)<<" euclid\n";
    return 0;*/
    std::cerr<<"Initialization done, running dim. reduction\n";
    if (nlrmode=="LLE") lleopts.mode=LLE;
    else if (nlrmode=="LLTE") lleopts.mode=LLTE;
    else if (nlrmode=="HLLE") lleopts.mode=HLLE;
    else if (nlrmode=="WLLE") { lleopts.mode=LLE; lleopts.nlopts.kw=nn; }
    else if (nlrmode=="WHLLE") { lleopts.mode=HLLE; lleopts.nlopts.kw=nn; }
    else if (nlrmode=="WLLTE") { lleopts.mode=LLTE; lleopts.nlopts.kw=nn; }
    else ERROR("Unsupported NLDR mode. Use LLE or LLTE.");
    
    NLDRLLE(mpoints,nlproj,lleopts,llereport);
    
    
    nlproj.get_points(hplist,lplist);
    
    if (fverb || fveryverb)
    {
        std::cout << " ######################## LLE REPORT ###################\n";
        std::cout << " # Error in fitting HD points: "<<llereport.hd_error<<"\n";
        std::cout << " # Small Eigenvalues of M: \n # ";
        for (int i=0; i<llereport.deval.size(); ++i) std::cout<<llereport.deval[i]<<" "; 
        if (fveryverb) std::cout<<"("<<llereport.dp1eval<<")"; std::cout <<"\n";
        std::cout << " # Error in fitting LD points: "<<llereport.ld_error<<"\n";
        std::cout << " # y1 .. yd "<<(fveryverb?" hd_error ld_error ":"")<<"\n";
    }
    
    for (int i=0; i<lplist.size(); i++)
    {
        for (int h=0; h<d; h++)  std::cout<<lplist[i][h]<<" ";
        if (fveryverb) std::cout<<llereport.hd_errors[i]<<" "<<llereport.ld_errors[i]<<" ";
        std::cout<<std::endl;
    }

    return 0;
}