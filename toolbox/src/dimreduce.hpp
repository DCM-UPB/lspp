#ifndef __DIMREDUCE_H
#define __DIMREDUCE_H 0
#include "tbdefs.hpp"
#include "matrix-full.hpp"
#include "matrix-crs.hpp"

namespace toolbox {
    
class NLDRLLEOptions; class NLDRLLEReport; 
class NLDRLLTEOptions; class NLDRLLTEReport; 
class NLDRWLLEOptions; class NLDRWLLEReport; 
class NLDRNeighborList;

class NLDROptions {
public:
    std::map<std::string,double> pars;
}; 

class NLDRProjection {
    friend void NLDRLLE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLEOptions& opts, NLDRLLEReport& report);
    
private:
    NLDROptions opts;
    unsigned long D, d, n;
    FMatrix<double> P, p;
    
public:
    void set_options(const NLDROptions& nopts)
    { opts=nopts; }
    void get_options(NLDROptions& nopts)
    { nopts=opts; }
    void set_points(const std::valarray<std::valarray<double> >& nP, const std::valarray<std::valarray<double> >& np)
    {
#ifdef DEBUG
        if ((n=nP.size())<1 || (D=nP[0].size())<1) ERROR("Hi-dimensional array has inconsistent sizes.");
        if (np.size()!=n || (d=np[0].size())<1 || d>D) ERROR("Low-dimensional array has inconsistent sizes.");
#endif
        P=nP; p=np; 
    }
    void get_points(std::valarray<std::valarray<double> >& nP, std::valarray<std::valarray<double> >& np) 
    { 
        nP.resize(n); np.resize(n); 
        for (unsigned long i=0; i<n; ++i) { nP[i].resize(D); for (unsigned long h=0; h<D; ++h) nP[i][h]=P(i,h); }
        for (unsigned long i=0; i<n; ++i) { np[i].resize(d); for (unsigned long h=0; h<d; ++h) np[i][h]=p(i,h); }
    }
};

/*! A collection of metric functions */
//null-metric class
class NLDRMetric {
friend class NLDRNeighborList;
private: 
    virtual double pdist(const double* a, const double* b, unsigned long n) const { return 0.0; };
    virtual void pdiff(const double* a, const double* b, double* c, unsigned long n) const { for (unsigned long i=0; i<n; ++i) c[i]=a[i]-b[i]; };
public:
    double dist(const std::valarray<double>& a, const std::valarray<double>& b) const
    { 
#ifdef DEBUG
        if (a.size()!=b.size()) ERROR("Vector size mismatch in distance.");
#endif
        return pdist(&a[0],&b[0],a.size()); 
    }
    void diff(const std::valarray<double>& a, const std::valarray<double>& b, std::valarray<double>& c) const
    {
#ifdef DEBUG
        if (a.size()!=b.size()) ERROR("Vector size mismatch in distance.");
#endif
        c.resize(a.size()); pdiff(&a[0], &b[0], &c[0], a.size()); 
    }
};

class NLDRMetricEuclid: public NLDRMetric {
private: 
    double pdist(const double* a, const double* b, unsigned long d) const;
};

class NLDRMetricPBC: public NLDRMetric {
private: 
    double pdist(const double* a, const double* b, unsigned long d) const;
    void pdiff(const double* a, const double* b, double* c, unsigned long n) const;
public:
    std::valarray<double> periods;
    NLDRMetricPBC() : periods() {}
    NLDRMetricPBC(const NLDRMetricPBC& no) : periods(no.periods) {}
    NLDRMetricPBC& operator= (const NLDRMetricPBC& no) 
    { if (&no==this) return *this; periods.resize(no.periods.size()); periods=no.periods; }
};

enum NLDRNeighborGreediness { NLDRGreedy, NLDRLiberal, NLDRAsym }; 
class NLDRNeighborOptions{
public:
    NLDRMetric *ometric;
    unsigned long kw;
    unsigned long maxneigh; double cutoff; NLDRNeighborGreediness greediness;
 
    NLDRNeighborOptions() : ometric(NULL), kw(0), maxneigh(4), greediness(NLDRGreedy), cutoff(0.0) {}
};

class NLDRNeighbor { public: unsigned long j; double d; NLDRNeighbor(unsigned long nj=0, double nd=0.0) :j(nj),d(nd) {}};
inline bool operator < (const NLDRNeighbor &lhs, const NLDRNeighbor &rhs) { return lhs.d< rhs.d; }

class NLDRNeighborList {
    friend void NLDRLLE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLEOptions& opts, NLDRLLEReport& report);
private:
    NLDRNeighborOptions opts;
    std::valarray<NLDRNeighbor> nlist;
    std::valarray<unsigned long>  npoint;
    void nlbuildup(const FMatrix<double>& points);
    NLDRNeighbor& rneigh(unsigned long i, unsigned long j); 
    
public:
    NLDRNeighborList& operator=(const NLDRNeighborList& nn)
    {
        if (&nn==this) return *this;
        nlist.resize(nn.nlist.size()); nlist=nn.nlist;
        npoint.resize(nn.npoint.size()); npoint=nn.npoint;
        opts=nn.opts;
    }
    NLDRNeighborList(const NLDRNeighborOptions& nopts=NLDRNeighborOptions()) : opts(nopts) {}
    NLDRNeighborList(const std::valarray<std::valarray<double> >& points, 
                     const NLDRNeighborOptions& nopts=NLDRNeighborOptions()) { Build(points, nopts); }
    NLDRNeighborList(const FMatrix<double>& points, 
                     const NLDRNeighborOptions& nopts=NLDRNeighborOptions()) : opts(nopts) { nlbuildup(points); }
    void Build(const FMatrix<double>& lpoints) { nlbuildup(lpoints); }
    void Build(const FMatrix<double>& lpoints, const NLDRNeighborOptions& nopts) { opts=nopts; nlbuildup(lpoints); }
    void RemoveLonesome(std::valarray<unsigned long>& llone);
    
    unsigned long size() const {return npoint.size()-1;}
    unsigned long nneigh(unsigned long i) const 
    {
#ifdef DEBUG
        if (i>=npoint.size()) ERROR("Index out of bounds in neighbor list");
#endif
        return npoint[i+1]-npoint[i]; 
    } 
    unsigned long index(unsigned long i, unsigned long j) const;
    double dist(unsigned long i, unsigned long j) const;  
    NLDRNeighbor neigh(unsigned long i, unsigned long j) const; 
};

enum NLDRLLEMode { LLE, LLTE, HLLE }; 
class NLDRLLEOptions
{
public:
    NLDRNeighborOptions nlopts; NLDRLLEMode mode; bool verbose; 
    unsigned long lowdim, dimts; double smooth; bool rmlonesome;
    NLDRLLEOptions() : nlopts(), mode(LLE), verbose(false), 
                   lowdim(2), dimts(2), smooth(-1e-5), rmlonesome(false) {}
};

class NLDRLLEReport
{
public:
    std::valarray<double> hd_errors; double hd_error;
    std::valarray<double> ld_errors; double ld_error;
    std::valarray<double> deval;  double dp1eval;

    NLDRLLEReport& operator=(const NLDRLLEReport& nr)
    {
        if (&nr==this) return *this;
        hd_errors.resize(nr.hd_errors.size()); hd_errors=nr.hd_errors;  hd_error=nr.hd_error;
        ld_errors.resize(nr.ld_errors.size()); ld_errors=nr.ld_errors;  ld_error=nr.ld_error;
        deval.resize(nr.deval.size()); deval=nr.deval;  dp1eval=nr.dp1eval;
        return *this;
    }
};

void NLDRLLE(FMatrix<double>& points, NLDRProjection& proj, const NLDRLLEOptions& opts, NLDRLLEReport& report);

}; //ends namespace toolbox
#endif //ends #ifndef __DIMREDUCE_H