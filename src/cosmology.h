#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include "cosmology.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spline.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cstdio>
#include <cstdlib>
#include <ctime>

struct cosmo
{
    double  Om0,Oml,Omb,hval,th,s8,nspec;
};

double findmvir(double, void*);
double findcmarch(double, void*);

class cosmology
{
    protected:
    /// Variables
    double Omega0,Omegal,Omegab,h,theta,sigma8,ns,d2norm,t0,rho_crit_0,facmtor;
    double fconc;

    bool verbose;

    private:
    /// Some constants
    double kmpspMpctoGyr; //=977.813952;
    double gee          ; //=4.2994e-9;
    double c            ; //=299792.458;
    double e            ; //=2.71828183;


    /// Options for various functions
    int opt_c;
    char cfile[500];
    //int opt_nm;
    //char nmfile[500];

    /// Private functions
    void initialize();
    double getcDel(double cvir, double z, double Delta);

    double c_Zhao(double M,double z);
    bool bool_zhao;
    double z_zhao;
    gsl_interp_accel *zhao_acc;
    gsl_spline *zhao_spline;
    void init_Zhao(double z);

    double c_gen(double M,double z);        // Concentration from file
    bool bool_gen;
    double z_gen;
    gsl_interp_accel *gen_acc;
    gsl_spline *gen_spline;
    void init_gen(double z);

    double xMmin,xMmax,cMmin,cMmax,dymin,dymax;

    bool bool_pe_rho_rdelta_phys_Zhao;
    gsl_interp_accel *pe_rho_rdelta_phys_Zhao_acc;
    gsl_spline *pe_rho_rdelta_phys_Zhao_spline;
    gsl_interp_accel *mah_Zhao_acc;
    gsl_spline *mah_Zhao_spline;
    void init_pe_rho_rdelta_phys_Zhao(double M,double z=0.0);
    double Mvir_for_pe,z_for_pe;

    // Auxiliary functions
    void getuuid(char *strUuid);
    
    
    //Friends
    friend double findmvir(double, void*);
    friend double findcmarch(double, void*);


    public:
    
        // Basic cosmology constructors
        cosmology();
        ~cosmology(); 
        cosmology(double Om0,double Oml,double Omb,double h,double thetaCMB,double sigma8,double nspec);
        cosmology(cosmo);
        void cosmo_free();

        double Delta_crit(double z);  //Delta_crit(z), Bryan and Norman '98
        double Omega(double z);       //Omega(z) 
        double Eofz(double z);        //Eofz(z) 


        //NFW profile related functions
        void modelNFWhalo_com_ext(double mDel,double z,double& Mvir,double& Rvir,double& cvir,double& Rdel,double& cDel, double Del);
        void getcDelta(double z,double cvir,double&cDel,double Del); 
        void modelNFWhalo_com(double mDel,double z,double& Mvir,double& Rvir,double& cvir,double Del);
        double rvir_from_mvir(double Mvir,double z); 
        void convertMdelta_Mdeltap(double z, double Mdelta, double Delta, double Deltap, double &Mdeltap, double &cdeltap);
        double getMvir(double M, double z, double Delta);
        double getM4rs(double M, double z);
        double dM4rs_dMphys(double mvir0,double z0,double z1,double z2);
        double dM4rs_dMvir(double mvir0,double z0,double z1,double z2);

        // Concentration parameter wrapper
        double conc(double M,double z);

        // Set concentration from file
        void set_cfile(char *cfile);
        void setfconc(double x);

        // Pseudo-evolution (forward/backward)
        void pevolve_fixed(double cdel,int opt,double z,double& cdelz,double& fdelz,double zstart=0.0);
        // Pseudo-evolution fraction
        double pe_fraction(double mvir0,double z0,double z1,double z2);

        // Output Mass accretion history from Zhao et al.
        void mah_Zhao(double M,double z=0.0);
        
};

#endif
