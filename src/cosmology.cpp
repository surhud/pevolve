///Cosmology routines
#include "cosmology.h"


/// march_params object
struct march_params
{
    double fac;
};


/// Object to pass cosmology, mDel, Delta and z
struct mvir_params
{
    cosmology *cptr;
    double *mDel;
    double *z;
    double *Delta;
};

/// Object to pass cosmology, cvir, Omega(z), Deltacrit(z)
struct cDel_params
{
    cosmology *cptr;
    double *cvir;
    double *omegaz;
    double *dcz;
    double *Delta;
};

/// mu(x) \int x/(1+x)^2 dx
double mu(double x){
    return log(1.+x)-x/(1.+x);
}

void cosmology::getuuid(char *strUuid){
    srand(time(NULL));
    sprintf(strUuid, "%x%x-%x-%x-%x-%x%x%x",
                rand(), rand(),                 // Generates a 64-bit Hex number
                rand(),                         // Generates a 32-bit Hex number
                ((rand() & 0x0fff) | 0x4000),   // Generates a 32-bit Hex number of the form 4xxx (4 indicates the UUID version)
                rand() % 0x3fff + 0x8000,       // Generates a 32-bit Hex number in the range [0x8000, 0xbfff]
                rand(), rand(), rand());        // Generates a 96-bit Hex number
    std::cout<<"# Temporary directory name to run mandc.x: "<<strUuid<<std::endl;
}

///Constructor WMAP3 cosmology
cosmology::cosmology()
{
    Omega0 = 0.27;
    Omegab = 0.047;
    Omegal = 1.-Omega0;
    h      = 0.7 ;
    theta  = 1.0093;///2.7;
    sigma8 = 0.82;
    ns     = 0.95;
    initialize();
    if(verbose) std::cout<<"# Cosmo constructor 1 called, initialized to Bolshoi\n";
}

/** Initialize the cosmology dependent routines
 *  and output to the log
 * */
void cosmology::initialize()
{

    // Some constants
    kmpspMpctoGyr=977.813952;
    gee          =4.2994e-9;
    c            =299792.458;
    e            =2.71828183;
    verbose      =true;

    /// Initialize options
    opt_c=1;

    bool_zhao=false;
    bool_gen=false;
    //bool_nm_gen=false;
    bool_pe_rho_rdelta_phys_Zhao=false;
    fconc=1.0;
    
    /// Fix some constants
    rho_crit_0=3.0e4/(8.*M_PI*gee);
    facmtor=(1./(4./3.*M_PI*Omega0*rho_crit_0));
    
    /// Output cosmology
    if(verbose){
    std::cout<<"# "<<"============================================================="<<std::endl;
    std::cout<<"# "<<"Cosmology initialised to: "<<std::endl;
    std::cout<<"# Om0="<<Omega0<<"\n" \
             <<"# Omb="<<Omegab<<"\n" \
             <<"# hub="<<h     <<"\n" \
             <<"# tht="<<theta <<"\n" \
             <<"# s8 ="<<sigma8<<"\n" \
             <<"# ns ="<<ns            \
             <<std::endl;
    std::cout<<"# t0 ="<<t0<<", ,rho_crit_0="<<rho_crit_0<<std::endl;
    std::cout<<"# "<<"============================================================="<<std::endl;
    }
    
}

///Constructor for cosmology by passing a cosmo struct
cosmology::cosmology(cosmo p)
{
    Omega0 = p.Om0;
    Omegab = p.Omb;
    Omegal = 1-p.Om0;
    h      = p.hval;
    theta  = p.th/2.7;
    sigma8 = p.s8;
    ns     = p.nspec;
    initialize();
    if(verbose) std::cout<<"# Cosmo constructor 2 called\n";
}

///Constructor for cosmology by passing all cosmological parameters explicitly
cosmology::cosmology(double om0, double oml, double omb, double xh ,double th, double s8, double nspec)
{
    Omega0 = om0;
    Omegab = omb;
    Omegal = oml;
    h      = xh;
    theta  = th/2.7;
    sigma8 = s8;
    ns     = nspec;
    initialize();
    if(verbose) std::cout<<"# Cosmo constructor 3 called\n";
}

///Destructor for cosmology object
cosmology::~cosmology()
{
    cosmo_free();
}


/// Expansion factor (z)
double cosmology::Eofz(double z)
{
    return sqrt(Omega0*pow(1.+z,3.0)+Omegal+(1-Omega0-Omegal)*pow(1.+z,2.0));
}

/// Omega(z)
double cosmology::Omega(double z)
{
    return Omega0*pow(1+z,3.)/pow(Eofz(z),2.);
}

/// Delta_crit(z) a'la Bryan and Norman '98
double cosmology::Delta_crit(double z)
{
    double result;

    double x=Omega(z)-1.;

    if (Omega0<1.&&Omegal==0.0)
    {
        result=18.0 * M_PI*M_PI + 60.0*x - 32.0 * x*x;
    } else if (fabs(Omega0+Omegal-1.0)<1.e-4)
    {
        result=18.0 * M_PI*M_PI + 82.0*x - 39.0 * x*x;
    } else
    {
        std::cout<<"# Aborting as Delta_Crit not defined for this cosmology."<<std::endl;
        exit(0);
    }

    return result;
}

void cosmology::cosmo_free(){
    if(verbose){
    std::cout<<"# Cosmo free destructor\n";
    }
    if(bool_zhao){
        gsl_interp_accel_free(zhao_acc);
        gsl_spline_free(zhao_spline);
        bool_zhao=false;
    }
    if(bool_gen){
        gsl_interp_accel_free(gen_acc);
        gsl_spline_free(gen_spline);
        bool_gen=false;
    }
}

void cosmology::set_cfile(char *file){
    strcpy(cfile,file);
    opt_c=2;
    fprintf(stdout,"# Using file %s for concentrations\n",cfile);
}

//// Routines which have to do with NFW Haloes

/// Concentration parameter wrapper
double cosmology::conc(double M, double z)
{
    double result;

    if(opt_c==1)
    {
        result=c_Zhao(M, z);
    } else if(opt_c==2)
    {
        result=c_gen(M, z);
    } else
    {
        fprintf(stderr,"Concentration parameter option not supported yet.\n");
        exit(0);
    }

    return result*fconc;
}

/// The concentration parameter form a file. Note that the mass should be Mvir in h^{-1} Msun
double cosmology::c_gen(double M, double z)
{
    /// First initialize if not initialized
    if (!bool_gen){
        init_gen(z);
    }else if(z!=z_gen){
        init_gen(z);
    }

    /// Spline interpolate
    double xM=log10(M);
    double concen;
    if(xM<xMmin){
        concen=cMmin+dymin*(xM-xMmin);
    }else{
        if(xM>xMmax){
            concen=cMmax+dymax*(xM-xMmax);
        }else{
            concen=gsl_spline_eval(gen_spline,xM,gen_acc);
        }
    }

    //std::cout<<xM<<" "<<concen<<std::endl;

    return pow(10.,concen);
}

void cosmology::init_gen(double z){

    /// Input file handler
    FILE *mcinp;

    /// Read the generic file and init spline
    mcinp=fopen(cfile,"r");
    double xx[200],yy[200];
    int fill=0;
    double xM, xc;
    //while(1){
    while ( fscanf(mcinp,"%lg %lg",&xM,&xc)==2 ){
	if(xc<0.0||fill>=200) break;
        xx[fill]=xM;
        yy[fill]=xc;
        fill++;
    }
    fclose(mcinp);
    gen_acc = gsl_interp_accel_alloc ();
    gen_spline = gsl_spline_alloc (gsl_interp_cspline, fill);
    gsl_spline_init (gen_spline,xx,yy,fill);

    /// Calculate the derivatives on the outer points
    xMmin=xx[0];
    xMmax=xx[fill-1];
    cMmin=yy[0];
    cMmax=yy[fill-1];
    dymin=(yy[1]-yy[0])/(xx[1]-xx[0]);
    dymax=(yy[fill-1]-yy[fill-2])/(xx[fill-1]-xx[fill-2]);

    /// Concentrations inited.
    if(verbose){
        std::cerr<<"# Generic concentration initialized for z: "<<z<<std::endl;
    }
    bool_gen=true;
    z_gen=z;
}

/// The concentration parameter ala Zhao et al. 2009. Note that the mass is Mvir in h^{-1} Msun
double cosmology::c_Zhao(double M, double z)
{
    /// First initialize if not initialized
    if (!bool_zhao){
        init_Zhao(z);
    }else if(z!=z_zhao){
        init_Zhao(z);
    }

    /// Spline interpolate
    double xM=log10(M);
    double concen;
    if(xM<xMmin){
        concen=cMmin+dymin*(xM-xMmin);
    }else{
        if(xM>xMmax){
            concen=cMmax+dymax*(xM-xMmax);
        }else{
            concen=gsl_spline_eval(zhao_spline,xM,zhao_acc);
        }
    }

    return pow(10.,concen);
}

/// This currently works only for Mvir, no enthusiasm of changing it for
/// Mdelta currently
void cosmology::mah_Zhao(double M,double z){
    FILE *mcinp;
    char dirname[100];
    getuuid(dirname);

    char *command;
    command =new char[2000];

    sprintf(command,"mkdir -p %s",dirname);
    system((const char *)command);
    
    char fname[105];
    sprintf(fname,"%s/tmpin",dirname);
    mcinp=fopen(fname,"w");

    fprintf(mcinp,"aumout\n");
    fprintf(mcinp,"%lg %lg \n",Omega0,Omegal);
    fprintf(mcinp,"1 \n");
    fprintf(mcinp,"%lg \n",h);
    fprintf(mcinp,"%lg \n",sigma8);
    fprintf(mcinp,"%lg \n",ns);
    fprintf(mcinp,"%lg %lg \n",Omegab,theta*2.7);
    fprintf(mcinp,"1 \n");
    fprintf(mcinp,"%lg \n",z);
    fprintf(mcinp,"%lg \n",log10(M));
    fclose(mcinp);

    sprintf(command,"cd %s; PATH=$PATH:.. mandc.x < tmpin >> zhao_log; rm s_aumout; rm sigma_aumout; mv mchistory* fileout;",dirname);
    system((const char*)command);

    delete [] command;

    /*
    system("PATH=$PATH:. mandc.x < tmpin >> tmp/zhao_log");
    system("rm s_aumout");
    system("rm sigma_aumout");
    system("mv mchistory* fileout");
    */

}

void cosmology::init_pe_rho_rdelta_phys_Zhao(double M,double z){

    FILE *mcinp;
    char dirname[100];
    getuuid(dirname);

    char *command;
    command =new char[2000];

    sprintf(command,"mkdir -p %s",dirname);
    system((const char *)command);
    
    char fname[105];
    sprintf(fname,"%s/tmpin",dirname);
    mcinp=fopen(fname,"w");
    //mcinp=fopen("tmpin","w");

    fprintf(mcinp,"aumout\n");
    fprintf(mcinp,"%lg %lg \n",Omega0,Omegal);
    fprintf(mcinp,"1 \n");
    fprintf(mcinp,"%lg \n",h);
    fprintf(mcinp,"%lg \n",sigma8);
    fprintf(mcinp,"%lg \n",ns);
    fprintf(mcinp,"%lg %lg \n",Omegab,theta*2.7);
    fprintf(mcinp,"1 \n");
    fprintf(mcinp,"%lg \n",z);
    fprintf(mcinp,"%lg \n",log10(M));
    fclose(mcinp);

    sprintf(command,"cd %s; PATH=$PATH:.. mandc.x < tmpin >> zhao_log; rm s_aumout; rm sigma_aumout; mv mchistory* fileout; awk '{if(NR!=1)print $1,$2,$3}' fileout > fileout.mvir.mah",dirname);
    system((const char *)command);

    /*
    system("PATH=$PATH:. mandc.x < tmpin >> tmp/zhao_log");
    system("rm s_aumout");
    system("rm sigma_aumout");
    system("mv mchistory* fileout");
    system("awk '{if(NR!=1)print $1,$2,$3}' fileout > fileout.mvir.mah");
    */

    
    /// What we want is a spline with rho_boundary(Rvir(z)) all physical 
    /// It allows us to calculate: \int_{R1}^{R2} 4 pi r^2 dr rho(zc,r)
    /// Rho(Rvir)=rhos/cvir/(1+cvir)^2
    /// rho(Rvir)=Mvir/(4 pi rs^3 mu(cvir) )/cvir/(1+cvir)^2

    double zred,mvir,cvir;
    double *xx,*yy,*zz,*xz,*cc;
    int fillmax=125;
    xx=new double [fillmax];
    yy=new double [fillmax];
    xz=new double [fillmax];
    zz=new double [fillmax];
    cc=new double [fillmax];

    sprintf(fname,"%s/fileout.mvir.mah",dirname);
    FILE *fp=fopen(fname,"r");
    int fill=0;
    for (int i=0;i<fillmax;i++){
        fscanf(fp,"%le %le %le",&zred,&mvir,&cvir);
	if(cvir<0.0 || fill>=fillmax){
	    break;
	}
        /// Obtain Rvir from mvir, convert to physical
        double rvir=rvir_from_mvir(mvir, zred);
        rvir=rvir/(1+zred);
        xx[fillmax-fill-1]=rvir;
        double rs=rvir/cvir;
        yy[fillmax-fill-1]=mvir/( 4*M_PI*pow(rs,3)*mu(cvir)*cvir*pow(1.+cvir,2) )*4*M_PI*rvir*rvir;
        zz[fill]=mvir;
        xz[fill]=zred;
        cc[fill]=cvir;
        fprintf(stdout,"%le %le %le\n",zred,mvir,cvir);
	fill++;
    }
    fclose(fp);

    // Shift the arrays yy and xx
    for(int i=0;i<fill;i++){
	yy[i]=yy[fillmax-fill+i];
	xx[i]=xx[fillmax-fill+i];
    }


    pe_rho_rdelta_phys_Zhao_acc=gsl_interp_accel_alloc ();
    pe_rho_rdelta_phys_Zhao_spline=gsl_spline_alloc(gsl_interp_cspline,fill);
    gsl_spline_init(pe_rho_rdelta_phys_Zhao_spline,xx,yy,fill);
    mah_Zhao_acc=gsl_interp_accel_alloc ();
    mah_Zhao_spline=gsl_spline_alloc(gsl_interp_cspline,fill);
    gsl_spline_init(mah_Zhao_spline,xz,zz,fill);
    cvir_mah_Zhao_acc=gsl_interp_accel_alloc ();
    cvir_mah_Zhao_spline=gsl_spline_alloc(gsl_interp_cspline,fill);
    gsl_spline_init(cvir_mah_Zhao_spline,xz,cc,fill);

    Mvir_for_pe=M;
    z_for_pe=z;
    bool_pe_rho_rdelta_phys_Zhao=true;

    delete []xx;
    delete []yy;
    delete []zz;

    sprintf(command,"rm -rf %s",dirname);
    system((const char *)command);

    delete [] command;

}

double cosmology::pe_fraction_fwd(double mvir0,double z0,double z1){
    if(!bool_pe_rho_rdelta_phys_Zhao || mvir0!=Mvir_for_pe || z0!=z_for_pe){
        std::cout<<"# Initializing physical density profile for "<<mvir0<<" at z="<<z0<<std::endl<<std::endl;
        init_pe_rho_rdelta_phys_Zhao(mvir0,z0);
    }
    if(z0>=z1){
        std::cout<<"# z0 should be smaller than z1 and z2\n";
        return 0;
    }
    /// Calculate the total mass evolution
    double mvir1=gsl_spline_eval(mah_Zhao_spline,z1,mah_Zhao_acc);

    double rvir1=rvir_from_mvir(mvir1,z1);
    rvir1=rvir1/(1+z1);
    double cvir1=gsl_spline_eval(cvir_mah_Zhao_spline,z1,cvir_mah_Zhao_acc);
    double rs1=rvir1/cvir1;

    double rvir0=rvir_from_mvir(mvir0,z0);
    rvir0=rvir0/(1+z0);

    /// Calculate the pseudo-evolution component
    double Mpe=mvir1/mu(cvir1)*( mu(rvir0/rs1) - mu(rvir1/rs1)  );

    return Mpe/(mvir0-mvir1);
}


double cosmology::pe_fraction(double mvir0,double z0,double z1,double z2){
    if(!bool_pe_rho_rdelta_phys_Zhao || mvir0!=Mvir_for_pe || z0!=z_for_pe){
        std::cout<<"# Initializing physical density profile for "<<mvir0<<" at z="<<z0<<std::endl<<std::endl;
        init_pe_rho_rdelta_phys_Zhao(mvir0,z0);
    }
    if(z0>z1 || z0>z2){
        std::cout<<"# z0 should be smaller than z1 and z2\n";
        return 0;
    }
    if(z2<z1){
        std::cout<<"# z1 should be smaller than z2\n";
        return 0;
    }
    /// Calculate the total mass evolution
    double Mvir1=gsl_spline_eval(mah_Zhao_spline,z1,mah_Zhao_acc);
    double Mvir2=gsl_spline_eval(mah_Zhao_spline,z2,mah_Zhao_acc);
    //mofz=(Mvir1+Mvir2)/2.;

    double Rvir1=rvir_from_mvir(Mvir1,z1);
    Rvir1=Rvir1/(1+z1);
    double Rvir2=rvir_from_mvir(Mvir2,z2);
    Rvir2=Rvir2/(1+z2);

    /// Calculate the pseudo-evolution component
    double Mpe=gsl_spline_eval_integ(pe_rho_rdelta_phys_Zhao_spline,Rvir2,Rvir1,pe_rho_rdelta_phys_Zhao_acc);
    //std::cout<<"DEBUG: "<<Mpe<<" "<<Mvir2-Mvir1<<std::endl;

    //printf("DEBUG : %e %e %e %e \n",Mvir1,Rvir1,Mvir2,Rvir2);
    //exit(101);

    return Mpe/(Mvir1-Mvir2);
}

double cosmology::dMcaustic_dMvir(double mvir0,double z0,double z1,double z2){
    if(!bool_pe_rho_rdelta_phys_Zhao || mvir0!=Mvir_for_pe || z0!=z_for_pe){
        std::cout<<"# Initializing physical density profile for "<<mvir0<<" at z="<<z0<<std::endl<<std::endl;
        init_pe_rho_rdelta_phys_Zhao(mvir0,z0);
    }
    if(z0>z1 || z0>z2){
        std::cout<<"# z0 should be smaller than z1 and z2\n";
        return 0;
    }
    if(z2<z1){
        std::cout<<"# z1 should be smaller than z2\n";
        return 0;
    }
    /// Calculate the total mass evolution
    double Mvir1=gsl_spline_eval(mah_Zhao_spline,z1,mah_Zhao_acc);
    double Mvir2=gsl_spline_eval(mah_Zhao_spline,z2,mah_Zhao_acc);
    //mofz=(Mvir1+Mvir2)/2.;

    /*
    double Rvir1=rvir_from_mvir(Mvir1,z1);
    Rvir1=Rvir1/(1+z1);
    double Rvir2=rvir_from_mvir(Mvir2,z2);
    Rvir2=Rvir2/(1+z2);
    */

    double dMvir=Mvir1-Mvir2;

    // I am assuming the same dlogMdloga for these objects, is this
    //reasonable?! Should not be that wrong an assumption
    double dlogMdloga=-(log10(Mvir1)-log10(Mvir2))/( log10(1+z1) - log10(1+z2) );

    // First calculate concentration of these halos
    double conc1=conc(Mvir1,z1);
    double conc2=conc(Mvir2,z2);

    // Get the c200_mean of these halos
    double c200m_1=getcDel(conc1,z1,200.0);
    double c200m_2=getcDel(conc2,z2,200.0);

    // Normalization factor for the Rt-R200m relation
    double norm_1=(1.+pow(Omega(z1)/Omega0,0.5))/2.;
    double norm_2=(1.+pow(Omega(z2)/Omega0,0.5))/2.;

    double c_caustic_1=c200m_1*norm_1*( 0.62+1.18*exp(-dlogMdloga/1.5) );
    double c_caustic_2=c200m_2*norm_2*( 0.62+1.18*exp(-dlogMdloga/1.5) );

    double Mcaustic1=Mvir1*mu(c_caustic_1)/mu(conc1);
    double Mcaustic2=Mvir2*mu(c_caustic_2)/mu(conc2);

    double dMcaustic=(Mcaustic1-Mcaustic2);
    //fprintf(stderr,"# Mvir1:%e Mvir2:%e Mpe:%e dMvir:%e M4rs1:%e M4rs2:%e dM4rs:%e\n",Mvir1,Mvir2,Mpe,dMvir,M4rs1,M4rs2,dM4rs);

    if(dMcaustic<0)
        return -1.0;
    else
        return dMcaustic/dMvir;
}


double cosmology::dM4rs_dMvir(double mvir0,double z0,double z1,double z2){
    if(!bool_pe_rho_rdelta_phys_Zhao || mvir0!=Mvir_for_pe || z0!=z_for_pe){
        std::cout<<"# Initializing physical density profile for "<<mvir0<<" at z="<<z0<<std::endl<<std::endl;
        init_pe_rho_rdelta_phys_Zhao(mvir0,z0);
    }
    if(z0>z1 || z0>z2){
        std::cout<<"# z0 should be smaller than z1 and z2\n";
        return 0;
    }
    if(z2<z1){
        std::cout<<"# z1 should be smaller than z2\n";
        return 0;
    }
    /// Calculate the total mass evolution
    double Mvir1=gsl_spline_eval(mah_Zhao_spline,z1,mah_Zhao_acc);
    double Mvir2=gsl_spline_eval(mah_Zhao_spline,z2,mah_Zhao_acc);
    //mofz=(Mvir1+Mvir2)/2.;

    double Rvir1=rvir_from_mvir(Mvir1,z1);
    Rvir1=Rvir1/(1+z1);
    double Rvir2=rvir_from_mvir(Mvir2,z2);
    Rvir2=Rvir2/(1+z2);

    /// Calculate the pseudo-evolution component
    //double Mpe=gsl_spline_eval_integ(pe_rho_rdelta_phys_Zhao_spline,Rvir2,Rvir1,pe_rho_rdelta_phys_Zhao_acc);
    //std::cout<<"DEBUG: "<<Mpe<<" "<<Mvir2-Mvir1<<std::endl;
    
    double dMvir=Mvir1-Mvir2;

    // First calculate concentration of these halos
    double conc1=conc(Mvir1,z1);
    double conc2=conc(Mvir2,z2);

    double M4rs1=getM4rs(Mvir1,conc1);
    double M4rs2=getM4rs(Mvir2,conc2);

    double dM4rs=(M4rs1-M4rs2);
    //fprintf(stderr,"# Mvir1:%e Mvir2:%e Mpe:%e dMvir:%e M4rs1:%e M4rs2:%e dM4rs:%e\n",Mvir1,Mvir2,Mpe,dMvir,M4rs1,M4rs2,dM4rs);

    if(dM4rs<0)
        return -1.0;
    else
        return dM4rs/dMvir;
}


double cosmology::dM4rs_dMphys(double mvir0,double z0,double z1,double z2){
    if(!bool_pe_rho_rdelta_phys_Zhao || mvir0!=Mvir_for_pe || z0!=z_for_pe){
        std::cout<<"# Initializing physical density profile for "<<mvir0<<" at z="<<z0<<std::endl<<std::endl;
        init_pe_rho_rdelta_phys_Zhao(mvir0,z0);
    }
    if(z0>z1 || z0>z2){
        std::cout<<"# z0 should be smaller than z1 and z2\n";
        return 0;
    }
    if(z2<z1){
        std::cout<<"# z1 should be smaller than z2\n";
        return 0;
    }
    /// Calculate the total mass evolution
    double Mvir1=gsl_spline_eval(mah_Zhao_spline,z1,mah_Zhao_acc);
    double Mvir2=gsl_spline_eval(mah_Zhao_spline,z2,mah_Zhao_acc);
    //mofz=(Mvir1+Mvir2)/2.;

    double Rvir1=rvir_from_mvir(Mvir1,z1);
    Rvir1=Rvir1/(1+z1);
    double Rvir2=rvir_from_mvir(Mvir2,z2);
    Rvir2=Rvir2/(1+z2);

    /// Calculate the pseudo-evolution component
    double Mpe=gsl_spline_eval_integ(pe_rho_rdelta_phys_Zhao_spline,Rvir2,Rvir1,pe_rho_rdelta_phys_Zhao_acc);
    //std::cout<<"DEBUG: "<<Mpe<<" "<<Mvir2-Mvir1<<std::endl;
    
    double dMphys=Mvir1-Mvir2-Mpe;

    // First calculate concentration of these halos
    double conc1=conc(Mvir1,z1);
    double conc2=conc(Mvir2,z2);

    double M4rs1=getM4rs(Mvir1,conc1);
    double M4rs2=getM4rs(Mvir2,conc2);

    double dM4rs=(M4rs1-M4rs2);
    //fprintf(stderr,"# Mvir1:%e Mvir2:%e Mpe:%e dMphys:%e M4rs1:%e M4rs2:%e dM4rs:%e\n",Mvir1,Mvir2,Mpe,dMphys,M4rs1,M4rs2,dM4rs);

    if(dM4rs<0)
        return -1.0;
    else
        return dM4rs/dMphys;
}

double cosmology::getM4rs(double M, double conc){
    return M*mu(4.0)/mu(conc);
}


void cosmology::init_Zhao(double z){
    /// Create temporary input file
    FILE *mcinp;
    char dirname[100];
    getuuid(dirname);

    char *command;
    command =new char[2000];

    sprintf(command,"mkdir -p %s",dirname);
    system((const char *)command);
    
    char fname[105];
    sprintf(fname,"%s/tmpin",dirname);
    mcinp=fopen(fname,"w");

    fprintf(mcinp,"aumout\n");
    fprintf(mcinp,"%lg %lg \n",Omega0,Omegal);
    fprintf(mcinp,"1 \n");
    fprintf(mcinp,"%lg \n",h);
    fprintf(mcinp,"%lg \n",sigma8);
    fprintf(mcinp,"%lg \n",ns);
    fprintf(mcinp,"%lg %lg \n",Omegab,theta*2.7);
    fprintf(mcinp,"2 \n");
    fprintf(mcinp,"%lg \n",z);

    fclose(mcinp);

    /// Run Zhao code
    //printf("#Running Zhao et al. 2009 code to init cosmology");
    sprintf(command,"cd %s; PATH=$PATH:.. mandc.x < tmpin >> zhao_log; rm s_aumout; rm sigma_aumout; mv mc_aumout* tmpout; awk '{if(NR>1)print log($2)/log(10.), log($3)/log(10.) }' tmpout > c_Zhao",dirname);
    system((const char *)command);

    /*
    system("mkdir -p tmp");
    system("PATH=$PATH:. mandc.x < tmpin >> tmp/zhao_log");
    system("rm s_aumout");
    system("rm sigma_aumout");
    system("mv mc_aumout* tmpout");
    //system("sh proc tmpout");
    system("awk '{if(NR>1)print log($2)/log(10.), log($3)/log(10.) }' tmpout > tmpout2");
    system("cat tmpin >> tmp/zhao_ops");
    system("echo '#####' >> tmp/zhao_ops");
    system("cat tmpout >> tmp/zhao_ops");
    system("echo '#####' >>  tmp/zhao_ops");
    system("cat tmpout2 >>  tmp/zhao_ops");
    system("echo '#####' >>  tmp/zhao_ops");
    system("mv tmpout2 c_Zhao");
    system("rm tmpin tmpout ");
    */

    /// Read Zhao output "c_Zhao" and init spline
    sprintf(fname,"%s/c_Zhao",dirname);
    //std::cout<<"Reading from file "<<fname<<std::endl;
    mcinp=fopen(fname,"r");
    double xx[200],yy[200];
    int fill=0;
    double xM, xc;
    while ( fscanf(mcinp,"%lg %lg",&xM,&xc)==2 ){
	if(xc<0.0||fill>=200) break;
        xx[fill]=xM;
        yy[fill]=xc;
        fill++;
    }
    fclose(mcinp);
    zhao_acc = gsl_interp_accel_alloc ();
    zhao_spline = gsl_spline_alloc (gsl_interp_cspline, fill);
    gsl_spline_init (zhao_spline,xx,yy,fill);

    /// Calculate the derivatives on the outer points
    xMmin=xx[0];
    xMmax=xx[fill-1];
    cMmin=yy[0];
    cMmax=yy[fill-1];
    dymin=(yy[1]-yy[0])/(xx[1]-xx[0]);
    dymax=(yy[fill-1]-yy[fill-2])/(xx[fill-1]-xx[fill-2]);

    /// Zhao et al. concentrations initialised
    if(verbose){
        std::cerr<<"# Zhao concentrations initialized for z: "<<z<<std::endl;
    }
    bool_zhao=true;
    z_zhao=z;

    sprintf(command,"rm -rf %s",dirname);
    system((const char *)command);

    delete [] command;
}

/// The root of this function gives the Mvir for a particular MDel
double findmvir(double mvir, void *params)
{
    mvir_params c1 = *(mvir_params *) params;
    cosmology *c2;
    double *mDel;
    double *z;
    double *Delta;
    c2=c1.cptr;
    mDel=c1.mDel;
    z=c1.z;
    Delta=c1.Delta;

    double cvir=(*c2).conc(mvir,*z);
    double cDel=(*c2).getcDel(cvir,*z,*Delta);

    double fcvir=log(1.+cvir) - cvir/(1.+cvir);
    double fcDel=log(1.+cDel) - cDel/(1.+cDel);

    return *mDel/mvir-fcDel/fcvir;

}

/// Get the Mvir for a given MDel
double cosmology::getMvir(double M, double z,double Delta)
{
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double m_lo = 0.1*M, m_hi = 5.0*M;

    gsl_function F;
    mvir_params p;
    p.cptr = this;
    p.mDel = &M;
    p.z=&z;
    p.Delta=&Delta;

    F.function = &findmvir;
    F.params = &p;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, m_lo, m_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        m_lo = gsl_root_fsolver_x_lower (s);
        m_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (m_lo, m_hi,0, 1.e-6);

        if (status == GSL_SUCCESS)
        {
            //std::cout<<"# "<<"zcollapse:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return res;

}

/// The root of this function gives the cDel for a particular cvir
double findcDel(double cDel, void *params)
{
    cDel_params c1 = *(cDel_params *) params;
    cosmology *c2;
    double *cvir;
    double *omegaz;
    double *dcz;
    double *Delta;
    c2=c1.cptr;
    cvir=c1.cvir;
    omegaz=c1.omegaz;
    dcz=c1.dcz;
    Delta=c1.Delta;

    double fcvir=log(1.+*cvir) - *cvir/(1.+*cvir);
    double fcDel=log(1.+cDel) - cDel/(1.+cDel);

    double res= (*dcz)/(*Delta*(*omegaz))*pow((*cvir)/cDel,3.0) - fcvir/fcDel;
    //printf("%lg %lg %lg %lg %lg \n",fcvir,fcDel,*cvir,cDel,res);
    //fflush(stdout);
    return res;

}

/// Get the cDel for a given cvir
double cosmology::getcDel(double cvir, double z, double Delta)
{
    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double c_lo = 0.01*cvir, c_hi = 10.0*cvir;

    gsl_function F;
    cDel_params p;
    p.cptr = this;
    p.cvir = &cvir;
    double omz=Omega(z);
    double dcz=Delta_crit(z);
    p.omegaz=&omz; 
    p.dcz=&dcz;
    p.Delta=&Delta;

    F.function = &findcDel;
    F.params = &p;
   
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, c_lo, c_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        c_lo = gsl_root_fsolver_x_lower (s);
        c_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (c_lo, c_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
            //std::cout<<"# "<<"zcollapse:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    return res;

}

double cosmology::rvir_from_mvir(double Mvir, double z)
{

    double Rvir=0.169*pow(Mvir/1.e12,1./3.);
    Rvir*=pow(Delta_crit(z)/178.0,-1.0/3.0);
    Rvir*=pow(Eofz(z),-2./3.);

    Rvir=Rvir*(1.+z);
    return Rvir;
}

/// Radii are in comoving units here
void cosmology::modelNFWhalo_com(double mDel,double z,double &Mvir,
double &Rvir, double &cvir, double Delta)
{
    Mvir=getMvir(mDel,z,Delta);

    Rvir=0.169*pow(Mvir/1.e12,1./3.);
    Rvir*=pow(Delta_crit(z)/178.0,-1.0/3.0);
    Rvir*=pow(Eofz(z),-2./3.);

    cvir=conc(Mvir,z);

    Rvir=Rvir*(1.+z);

}

/// Radii are in comoving units here
void cosmology::getcDelta(double z, double cvir, double &cDel,double Del)
{

    double Delta=Del;

    cDel=getcDel(cvir,z,Delta);

}

/// Radii are in comoving units here
void cosmology::modelNFWhalo_com_ext(double mDel,double z,double &Mvir,
double &Rvir, double &cvir, double &RDel, double &cDel,double Del)
{
    //std::cout<<Delta_crit(z)<<std::endl;exit(0);
    Mvir=getMvir(mDel,z,Del);

    Rvir=0.169*pow(Mvir/1.e12,1./3.);
    Rvir*=pow(Delta_crit(z)/178.0,-1.0/3.0);
    Rvir*=pow(Eofz(z),-2./3.);

    double Delta=Del;
    RDel=pow( 1.0e12/(4./3.*Delta*3e4*0.3/(8*gee)),1./3.);
    RDel*=pow(mDel/1.e12,1./3.);
    RDel*=pow(Omega(z)/0.3,-1.0/3.0);
    RDel*=pow(Eofz(z),-2./3.);

    cvir=conc(Mvir,z);

    cDel=getcDel(cvir,z,Delta);

    Rvir=Rvir*(1.+z);
    RDel=RDel*(1.+z);

}

/// The root of this function gives the concentration for a halo
/// evolving in isolation
double findcmarch(double cdelz, void *params)
{
    march_params c1 = *(march_params *) params;
    double fac;
    fac=c1.fac;

    double res= fac - pow(cdelz,3)/mu(cdelz);
    //std::cerr<<cdelz<<" "<<res<<std::endl;
    return res;

}

void cosmology::pevolve_fixed(double cdel, int opt, double z, double
&cdelz, double &fdelz, double zstart){

    double fac;
    if(opt==1){
        fac=pow(cdel*(1.+zstart)/(1.+z),3)/mu(cdel);
    }else if(opt==2){
        fac=pow(cdel,3)/mu(cdel);
        fac=fac*pow(Eofz(zstart)/Eofz(z),2);
        fac=fac*Delta_crit(zstart)/Delta_crit(z);
    }else if(opt==3){
        fac=pow(cdel,3)/mu(cdel);
        fac=fac*pow(Eofz(zstart)/Eofz(z),2);
    }else{
        fprintf(stderr,"Option %d not supported yet, bailing out...",opt);
        exit(100);
    }

    int status;
    int iter = 0, max_iter = 100;
    double res;

    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    double c_lo = 0.01*cdel, c_hi = 100000.0*cdel;

    gsl_function F;
    march_params p;
    p.fac = fac;

    F.function = &findcmarch;
    F.params = &p;
   
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, c_lo, c_hi);

    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        res = gsl_root_fsolver_root (s);
        c_lo = gsl_root_fsolver_x_lower (s);
        c_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (c_lo, c_hi,0, 1e-6);

        if (status == GSL_SUCCESS)
        {
            //std::cout<<"# "<<"zcollapse:Brent converged after "<< iter<<" iterations"<<std::endl;
        }


    }while (status == GSL_CONTINUE && iter < max_iter);

    gsl_root_fsolver_free (s);

    cdelz=res;
    fdelz=mu(cdelz)/mu(cdel);

}

void cosmology::setfconc(double x){
    fconc=x;
}

void cosmology::convertMdelta_Mdeltap(double z, double Mdelta, double Delta, double Deltap, double &Mdeltap, double &cdeltap){
    double Mvir,Rvir,cvir;
    modelNFWhalo_com(Mdelta,z,Mvir,Rvir,cvir,Delta);

    cdeltap=getcDel(cvir,z,Deltap);
    Mdeltap=Mvir*mu(cdeltap)/mu(cvir);
    
}


/*
/// The mass function from a generic file. Note that the mass should
/// corespond to the same definition as desired
double cosmology::nofm_gen(double M, double z)
{
    /// First initialize if not initialized
    if (!bool_nm_gen){
        init_nm_gen(z);
    }else if(z!=z_nm_gen){
        init_nm_gen(z);
    }

    /// Spline interpolate
    double xM=log10(M);
    double nofm;
    int option=-1;
    if(xM<xMmin_nm){
        nofm=nMmin+dymin_nm*(xM-xMmin_nm);
        option=1;
    }else{
        if(xM>xMmax_nm){
            nofm=nMmax+dymax_nm*(xM-xMmax_nm);
        option=2;
        }else{
            nofm=gsl_spline_eval(gen_nm_spline,xM,gen_nm_acc);
        option=3;
        }
    }

    //std::cout<<xM<<" "<<nofm<<" "<<option<<std::endl;

    return pow(10.,nofm);
}

void cosmology::init_nm_gen(double z){

    /// Input file handler
    FILE *nminp;

    /// Read the generic file and init spline
    nminp=fopen(nmfile,"r");
    double xx[200],yy[200];
    int fill=0;
    double xM, xnofm;
    //while(1){
    while ( fscanf(nminp,"%lg %lg",&xM,&xnofm)==2 ){
        xx[fill]=log10(xM);
        yy[fill]=log10(xnofm);
        fill++;
    }
    fclose(nminp);
    gen_nm_acc = gsl_interp_accel_alloc ();
    gen_nm_spline = gsl_spline_alloc (gsl_interp_cspline, fill);
    gsl_spline_init (gen_nm_spline,xx,yy,fill);

    /// Calculate the derivatives on the outer points
    xMmin_nm=xx[0];
    xMmax_nm=xx[fill-1];
    nMmin=yy[0];
    nMmax=yy[fill-1];
    dymin_nm=(yy[1]-yy[0])/(xx[1]-xx[0]);
    dymax_nm=(yy[fill-1]-yy[fill-2])/(xx[fill-1]-xx[fill-2]);

    /// Mass function inited.
    if(verbose){
        std::cerr<<"# Generic Mass function initialized for z: "<<z<<std::endl;
    }
    bool_nm_gen=true;
    z_nm_gen=z;
}
*/
