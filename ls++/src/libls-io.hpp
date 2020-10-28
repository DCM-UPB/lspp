//I/O for LSHelper class
#ifdef __LIBLS_HPP
namespace toolbox{
    std::istream& operator>>(std::istream& istr, LSFMode& lsf)
    {
        std::string sn;
        istr >> sn; if (!istr.good()) ERROR("Error while reading from stream.");
        if (sn=="fastpoly") lsf=lsf_fastpoly;
        else if (sn=="niextrapolation") lsf=lsf_niextrapolation;
        else if (sn=="niwarmstart") lsf=lsf_niwarmstart;
        else if (sn=="autoqbar") lsf=lsf_autoqbar;
        else if (sn=="optk") lsf=lsf_optk;
        else if (sn=="autohamspec") lsf=lsf_autohamspec;
        else if (sn=="niextra2") lsf=lsf_niextra2;
        else if (sn=="niextra3") lsf=lsf_niextra3;
        else if (sn=="autolength") lsf=lsf_autolength;
        else if (sn=="keepinverse") lsf=lsf_keepinverse;
        else ERROR("Invalid flag for LSFMode enum.");
        return istr;
    }

    std::ostream& operator<<(std::ostream& ostr, const LSFMode& lsf)
    {
        switch (lsf) {
            case lsf_fastpoly: ostr<<"fastpoly ";  break;
            case lsf_niextrapolation: ostr<<"niextrapolation ";  break;
            case lsf_niwarmstart: ostr<<"niwarmstart ";  break;
            case lsf_autoqbar: ostr<<"autoqbar ";  break;
            case lsf_optk: ostr<<"optk ";  break;
            case lsf_autohamspec: ostr<<"autohamspec ";  break;
            case lsf_niextra2: ostr<<"niextra2 ";  break;
            case lsf_niextra3: ostr<<"niextra3 ";  break;
            case lsf_autolength: ostr<<"autolength ";  break;
            case lsf_keepinverse: ostr<<"keepinverse ";  break;
        }
        return ostr;
    };

    template<> bool IField<LSOptions>::operator>>(std::ostream& ostr) const {
        IOMap iom; std::vector<LSFMode > lsfm;
        iom.insert(value.beta,"beta");
        iom.insert(value.mu,"mu");
        iom.insert(value.hs_min,"ham_min");
        iom.insert(value.hs_max,"ham_max");
        iom.insert(value.accu,"accuracy");
        iom.insert(value.mops,"matrix_options");
        lsfm.resize(0); 
    
        if (value.flags & lsf_fastpoly) lsfm.push_back(lsf_fastpoly);
        if (value.flags & lsf_niextrapolation) lsfm.push_back(lsf_niextrapolation);
        if (value.flags & lsf_niwarmstart) lsfm.push_back(lsf_niwarmstart);
        if (value.flags & lsf_autoqbar) lsfm.push_back(lsf_autoqbar);
        if (value.flags & lsf_optk) lsfm.push_back(lsf_optk);
        if (value.flags & lsf_autohamspec) lsfm.push_back(lsf_autohamspec);
        if (value.flags & lsf_autohamspec) lsfm.push_back(lsf_autohamspec);
        if (value.flags & lsf_keepinverse) lsfm.push_back(lsf_keepinverse);
        iom.insert(lsfm,"flags");
    
        iom.insert(value.qbar,"q_bar");
        iom.insert(value.p,"p");
        iom.insert(value.ham_file,"ham_file");
        iom.insert(value.rho_file,"rho_file");
        iom.insert(value.frestart,"restart");
        ostr<<name<<" {\n"<<iom<<"\n}\n";
        return true;
    }

    template<> bool IField<const LSOptions>::operator>>(std::ostream& ostr) const {
        IOMap iom; std::vector<LSFMode > lsfm;
        iom.insert(value.beta,"beta");
        iom.insert(value.mu,"mu");
        iom.insert(value.hs_min,"ham_min");
        iom.insert(value.hs_max,"ham_max");
        iom.insert(value.accu,"accuracy");
        iom.insert(value.mops,"matrix_options");

        lsfm.resize(0);
 
        if (value.flags & lsf_fastpoly) lsfm.push_back(lsf_fastpoly);
        if (value.flags & lsf_niextrapolation) lsfm.push_back(lsf_niextrapolation);
        if (value.flags & lsf_niwarmstart) lsfm.push_back(lsf_niwarmstart);
        if (value.flags & lsf_autoqbar) lsfm.push_back(lsf_autoqbar);
        if (value.flags & lsf_optk) lsfm.push_back(lsf_optk);
        if (value.flags & lsf_autohamspec) lsfm.push_back(lsf_autohamspec);
        if (value.flags & lsf_autohamspec) lsfm.push_back(lsf_autohamspec);
        if (value.flags & lsf_keepinverse) lsfm.push_back(lsf_keepinverse);
        iom.insert(lsfm,"flags");
        iom.insert(value.qbar,"q_bar");
        iom.insert(value.p,"p");
        iom.insert(value.ham_file,"ham_file");
        iom.insert(value.rho_file,"rho_file");
        iom.insert(value.frestart,"restart");
        ostr<<name<<" {\n"<<iom<<"\n}\n";
        return true;
    }

    template<> bool IField<LSOptions>::operator<<(std::istream& istr) {
        IOMap iom; std::vector<LSFMode > lsfm(0);
        iom.insert(value.beta,"beta",1.);
        iom.insert(value.mu,"mu",0.);
        iom.insert(value.hs_min,"ham_min",0.); 
        iom.insert(value.hs_max,"ham_max",0.);
        iom.insert(value.accu,"accuracy",1e-4);
        iom.insert(lsfm,"flags");
        iom.insert(value.qbar,"q_bar",0x1ul);
        iom.insert(value.p,"p",0x10000ul);
        iom.insert(value.ham_file,"ham_file",std::string(""));
        iom.insert(value.rho_file,"rho_file",std::string(""));
        iom.insert(value.frestart,"restart",false);
        iom.insert(value.mops,"matrix_options");
        iom << istr;
        value.flags=0; for (unsigned long i=0; i<lsfm.size(); ++i) value.flags |= lsfm[i]; 
        flags|=iff_set; flags&=~iff_nonvalid;
        IOMap::iterator it=iom.begin();
        for (;it!=iom.end();++it) 
        {
            if (!((it->second->flags & iff_optional) || (it->second->flags & iff_set))) { flags&=~iff_set; return false;}
            if (it->second->flags & iff_nonvalid) { flags|=iff_nonvalid; return false;}
        }
        return true;
    }
   class Test_Complex : public IFBase
    {
   private: 
    C2Matrix value;


   public: 
    Test_Complex(CMatrix var, const std::string& nname) : IFBase(nname), value(var) { flags |=iff_optional; };
    void MakeMat(CMatrix AH, int sz)
    {
     for (int i=0; i<sz; ++i)
     for (int j=0; j<sz; ++j)
     {
//        if (i = j) 
        AH(i,j)= -0.1;
     }
        value = AH; 
    }
   void print() 
    {
      std::cout << "Matrix_out --> Start!" << std::endl;
      std::cout << value << std::endl;
      std::cout << "Matrix_out --> End!" << std::endl;
    } 
 //   ~Test_Complex(); 
 
};
 

    template<> bool IField<LSHelper>::operator>>(std::ostream& ostr) const 
    {
        IOMap iom;
        iom.insert(value.pops,"options");
        if (value.pops.frestart==true) iom.insert(value.HAM,"HAM");
        if (value.pops.frestart==true) iom.insert(value.RHO,"RHO");
        ostr<<name<<" {\n"<<iom<<"\n}\n";
        return true;
    }

    template<> bool IField<const LSHelper>::operator>>(std::ostream& ostr) const 
    {
        IOMap iom;
        iom.insert(value.pops,"options");
        if (value.pops.frestart==true) iom.insert(value.HAM,"HAM");
        if (value.pops.frestart==true) iom.insert(value.RHO,"RHO");
        ostr<<name<<" {\n"<<iom<<"\n}\n";
        return true;
    }

    template<> bool IField<LSHelper>::operator<<(std::istream& istr)
    {
    //the input must have the options BEFORE the ham. data, if restart is required!
        int sz=5; 
        CMatrix AH(sz,sz);

        std::string label;
        while ( (istr>>label) && label!="options" );
    
        if (label!="options") ERROR("Input file must contain 'options' field first.");
        IField<LSOptions> lso(value.pops, "options");
        lso<<istr;

        if(lso.flags & iff_nonvalid) ERROR("Invalid options in input file.");
        if(!((lso.flags & iff_set)|(lso.flags & iff_optional))) ERROR("Required options missing in input file.");
        if (value.pops.frestart==true)
        {
        //we read ham & rho from this very same input
            IOMap iom;
            iom.insert(value.HAM,"HAM");
            iom.insert(value.RHO,"RHO");
        }
        else
        {
            std::ifstream ifs;
        //we look for other files
            if (value.pops.ham_file!="")
            {
                ifs.open(value.pops.ham_file.c_str());
                if (!ifs.good()) ERROR("Error opening hamiltonian file.");
                //std::cout << "Before IField!" << std::endl;
                //IField<CMatrix> ifham(value.HAM,"");
                //ifs >> ifham;
                Test_Complex ifham(value.HAM,"");
                ifham.MakeMat(AH,sz);
                ifham.print();
                //std::cout << "After IField!" << std::endl;
                if (ifham.flags & iff_nonvalid) ERROR("Error reading hamiltonian from file.");
                ifs.close();
            }
            else ERROR("No hamiltonian given, either as restart or as external file.");
        }

        //std::cout << "After Constr." << std::endl;
        std::ofstream ofs;
        ofs.open("reproduction_h");
        IField<CMatrix> iflh(value.HAM,"");
        ofs << iflh;
        ofs.close(); 
    
        value.init_all();
        return istr;
    }

} //ENDS NAMESPACE TOOLBOX
#endif

#ifdef __LIBLS_CHEB_HPP
namespace toolbox {
    template<> bool IField<LSChebHelper>::operator>>(std::ostream& ostr) const 
    {
        IOMap iom;
        iom.insert(value.pops,"options");
        if (value.pops.frestart==true) iom.insert(value.HAM,"HAM");
        if (value.pops.frestart==true) iom.insert(value.RHO,"RHO");
        ostr<<name<<" {\n"<<iom<<"\n}\n";
        return true;
    }

    template<> bool IField<const LSChebHelper>::operator>>(std::ostream& ostr) const 
    {
        IOMap iom;
        iom.insert(value.pops,"options");
        if (value.pops.frestart==true) iom.insert(value.HAM,"HAM");
        if (value.pops.frestart==true) iom.insert(value.RHO,"RHO");
        ostr<<name<<" {\n"<<iom<<"\n}\n";
        return true;
    }

    template<> bool IField<LSChebHelper>::operator<<(std::istream& istr)
    {
    //the input must have the options BEFORE the ham. data, if restart is required!
        bool fnoham=false;
        std::string label;
        while ( (istr>>label) && label!="options" );
    
        if (label!="options") ERROR("Input file must contain 'options' field first.");
        IField<LSOptions> lso(value.pops, "options");
        lso<<istr;

        if(lso.flags & iff_nonvalid) ERROR("Invalid options in input file.");
        if(!((lso.flags & iff_set)|(lso.flags & iff_optional))) ERROR("Required options missing in input file.");
        if (value.pops.frestart==true)
        {
        //we read ham & rho from this very same input
            IOMap iom;
            iom.insert(value.HAM,"HAM");
            iom.insert(value.RHO,"RHO");
        }
        else 
        {
            std::ifstream ifs;
        //we look for other files
            if (value.pops.ham_file!="")
            {
                std::pout<<"Reading hamiltonian from file: "<<value.pops.ham_file.c_str()<<"\n";
                ifs.open(value.pops.ham_file.c_str());
                if (!ifs.good()) ERROR("Error opening hamiltonian file.");
                IField<RMatrix> ifham(value.HAM,"");
                ifs >> ifham;
                if (ifham.flags & iff_nonvalid) ERROR("Error reading hamiltonian from file.");
                ifs.close();
            }
            else fnoham=true;
        }
    
        if (!fnoham)
        {
            value.HAM.setops(value.pops.mops);
            std::pout<<" # HAMILTONIAN DATA: trace: "<<trace(value.HAM)<<"  norminf: "<<norminf(value.HAM)<<"  normfrob: "<<normfrob(value.HAM)<<"\n";
//             std::cout << " ham p r c: " << value.HAM.wr << " " << value.HAM.wc << std::endl;
//
             std::ofstream ofs;
             ofs.open("reproduction_h");
             IField<RMatrix> iflh(value.HAM,"");
             ofs << iflh;
             ofs.close();

            value.init_all();
        }
        else WARNING("No hamiltonian data given in input file.");
        return istr;
    }

} //ENDS NAMESPACE TOOLBOX
#endif
