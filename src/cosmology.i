%module cosmology
%include typemaps.i
%apply int &OUTPUT { int & fA };
%apply double &OUTPUT { double & Mvir, double &Rvir, double& cvir,double&Rdel,double& cDel };
%apply double &OUTPUT { double & Mvir, double &Rvir, double& cvir };
%apply double &OUTPUT { double &cDel };
%apply double &OUTPUT { double& cdelz, double &fdelz };
%apply double &OUTPUT { double &mofz };
%apply double &OUTPUT { double &Mdeltap, double &cdeltap };
%feature("autodoc", 1);
%{
    #define SWIG_FILE_WITH_INIT
    #include "cosmology.h"
%}

%feature("docstring") cosmology::cosmology 
"Initializes cosmology object

:Parameters:

-   Omega0 : Matter density parameter
-   OmegaLambda : Dark energy density parameter
-   Omegab : Baryon density parameter
-   h : Hubble parameter
-   ThetaCMB : CMB temperature
-   sigma8 : sigma8
-   ns : power spectrum index

:Returns:

-   Cosmology object

Without any inputs, initializes to Bolshoi cosmology.

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> help(a)
"

%feature("docstring") cosmology::Eofz
"Returns the cosmological expansion function E(z)

:Parameters:

-   z : Redshift

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.Eofz(0.5)
    1.28111279753
"

%feature("docstring") cosmology::Omega 
"Returns the matter density parameter at redshift z

:Parameters:

-   z : Redshift

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.Omega(0.5)
    0.555217060168
"

%feature("docstring") cosmology::Delta_crit
"Returns the critical overdensity for collapse at redshift z
with respect to the critical density of the Universe a'la
Bryan and Norman 98

:Parameters:

-   z : Redshift

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.Delta_crit(0.0)
    97.0097792196
"

%feature("docstring") cosmology::set_cfile 
"Sets the file to be used for the concentration mass relation

:Parameters:

-   cfile : File name that contains cvir(Mvir) that should be read. The columns should be Mvir, cvir

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.set_cfile('myconcentration.dat')
"

%feature("docstring") cosmology::conc
"Returns concentration as a function of Mvir and z

:Parameters:

-  Mvir : Virial mass in hinv Msun
-  z : Redshift

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.conc(1.0E13,0.0)
    7.57896671931
"

%feature("docstring") cosmology::mah_Zhao
"Output the mass accretion history of a halo of mass M at
redshift z, output stored in file fileout

:Parameters:

-   M : Virial mass in hinv Msun
-   z : Redshift

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.mah_Zhao(1.E11,0.5)
"

%feature("docstring") cosmology::pe_fraction_fwd
"Output the pseudo-evolution fraction of halo of mass mvir0 at
redshift z0 between redshifts z0 and z1, such that z0<z1, assuming that the
density profile at z1 is unchanged

:Parameters:

-   mvir0 : Virial mass in hinv Msun
-   z0 : Redshift at which the mvir0 is identified
-   z1 : Initial redshift
    
:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.pe_fraction(1.E11,0.0,1.0)
"

%feature("docstring") cosmology::pe_fraction 
"Output the pseudo-evolution fraction of halo of mass Mvir0 at
redshift z0 between redshifts z1 and z2, z0<z1<z2

:Parameters:

-   mvir0 : Virial mass in hinv Msun
-   z0 : Redshift at which the mvir0 is identified
-   z1 : Initial redshift
-   z2 : Final redshift
    
:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.pe_fraction(1.E11,0.0,1.0,2.0)
    0.614891598653
"

%feature("docstring") cosmology::getMvir 
"Output the virial mass of a halo of mass M at redshift z but
defined as Delta with respect to the background density

:Parameters:

-   M : Mass (hinv Msun) defined as Delta times overdense with respect to background density
-   z : Redshift at which M is identified
-   Delta : The overdensity with respect to the background

:Returns:

-   Mvir : Virial mass in hinv Msun

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.getMvir(1.E11,0.0,200.0)
    88953346252.6
"

%feature("docstring") cosmology::getcDelta 
"Output the concentration, cDelta, of a halo defined as Delta
with respect to the background density at redshift z and which
has virial concentration equal to cvir.

:Parameters:

-   z : Redshift
-   cvir : The virial concentration
-   Del : The overdensity with respect to the background

:Returns:

-   cDelta : The concentration of this halo when defined as Del times overdense with respect to the background

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.getcDelta(0.0,10.0,200.)
    12.6784160959
"

%feature("docstring") cosmology::rvir_from_mvir
"Output the comoving virial radius of a halo of mass Mvir at
redshift z

:Parameters:

-   Mvir : Virial mass in hinv Msun
-   z : Redshift

:Returns:

-   Rvir : The comoving virial radius in hinv Mpc

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.rvir_from_mvir(1.E12,0.0)
    0.20689732071
"

%feature("docstring") cosmology::modelNFWhalo_com
"Output the virial mass, comoving virial radius, and virial
concentration of a halo of mass Mdelta at redshift z, defined
as Delta with respect to the background density

:Parameters:

-   mDel : Mass (hinv Msun) defined as Delta times overdense with respect to background density
-   z : Redshift
-   Delta : The overdensity

:Returns:

-   Mvir : The virial mass (hinv Msun)
-   Rvir : The comoving virial radius (hinv Mpc)
-   cvir : The virial concentration

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.modelNFWhalo_com(1.E12,0.0,200.0)
    [880330641353.4209, 0.198291214519584, 9.765326694230426]
"

%feature("docstring") cosmology::modelNFWhalo_com_ext
"Output the virial mass, comoving virial radius, virial
concentration, boundary radius Rdelta, and concentration
cDelta of a halo of mass Mdelta at redshift z, defined as
Delta with respect to the background density

:Parameters:

-   mDel : Mass (hinv Msun) defined as Delta times overdense with respect to background density
-   z : Redshift
-   Delta : The overdensity

:Returns:

-   Mvir : The virial mass (hinv Msun)
-   Rvir : The comoving virial radius (hinv Mpc)
-   cvir : The virial concentration
-   Rdel : The comoving boundary of halo Delta times overdense with respect to background density (hinv Mpc)
-   cdel : The concentration of halo Delta times overdense with respect to background density

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.modelNFWhalo_com(1.E12,0.0,200.0)
    [880330641353.4209, 0.198291214519584, 9.765326694230426, 0.2515830411841534, 12.386409188477131]
"

%feature("docstring") cosmology::pevolve_fixed
"Pseudo-evolution estimate for the mass (backward or forward):
Assume that the physical density profile of a peak with
concentration cdel at redshift zstart defined to be of type
opt 
remains fixed. Calculate its concentration at redshift z, and
the ratio of its mass to the mass at redshift zstart.

:Parameters:

-   cdel : concentration of halo 
-   opt : opt=1: Defined with respect to background density, opt=2: Defined to be virial mass, opt=3: Defined with respect to critical density
-   z: Redshift
-   zstart : The reference redshift at which the halo density profile is fixed

:Returns:

-   cdelz : The concentration of the halo at redshift z
-   fdelz : The ratio of the mass at redshift z to the mass at redshift zstart

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.pevolve_fixed(12.0,1,1.0,0.0)
    [5.019071157921484, 0.585350923353302]
"

%feature("docstring") cosmology::setfconc
"Multiply all concentrations by a factor x from here on

:Parameters:

-   x : multiplicative factor to the concentration of halo 
"

%feature("docstring") cosmology::convertMdelta_Mdeltap 
"Output the mass and concentration of a halo defined to be
Deltap times overdense with respect to background density if
it has mass Mdelta when defined as Delta times overdense with
respect to the background at redshift z

:Parameters:

-   z : Redshift
-   Mdelta : Mass (hinv Msun) defined as Delta times overdense with respect to background density
-   Delta : The overdensity
-   Deltap : The new overdensity

:Returns:

-   Mdeltap : The mass (hinv Msun) defined as Deltap times overdense with respect to background density
-   cdeltap : The concentration of halo defined as Deltap times overdense with respect to background density

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.convertMdelta_Mdeltap(0.0,1.0E12,200.,300.)
    [916906871432.2742, 10.512114412989405]
"

%feature("docstring") cosmology::dM4rs_dMphys 
"Output the ratio of the growth in M4rs and Mvir*(1-pe_fraction) between
redshifts z1 and z2, z0<z1<z2 for a halo of mass mvir0 at redshift z0 

:Parameters:

-   mvir0 : Virial mass in hinv Msun
-   z0 : Redshift at which the mvir0 is identified
-   z1 : Initial redshift
-   z2 : Final redshift
    
:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.dM4rs_dMphys(1.E11,0.0,1.0,2.0)
"

%feature("docstring") cosmology::Mcaustic_from_Mvir 
"Output Mcaustic for a given Mvir at redshift z0

:Parameters:

-   mvir0 : Virial mass in hinv Msun
-   z0 : Redshift at which the mvir0 is identified
    
:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.Mcaustic_from_Mvir(1.E11,0.0)
"

%feature("docstring") cosmology::dMcaustic_dMvir 
"Output the ratio of the growth in Mcaustic and Mvir between
redshifts z1 and z2, z0<z1<z2 for a halo of mass mvir0 at redshift z0 

:Parameters:

-   mvir0 : Virial mass in hinv Msun
-   z0 : Redshift at which the mvir0 is identified
-   z1 : Initial redshift
-   z2 : Final redshift
    
:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.dMcaustic_dMvir(1.E11,0.0,1.0,2.0)
"

%feature("docstring") cosmology::dM4rs_dMvir 
"Output the ratio of the growth in M4rs and Mvir between
redshifts z1 and z2, z0<z1<z2 for a halo of mass mvir0 at redshift z0 

:Parameters:

-   mvir0 : Virial mass in hinv Msun
-   z0 : Redshift at which the mvir0 is identified
-   z1 : Initial redshift
-   z2 : Final redshift
    
:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.dM4rs_dMvir(1.E11,0.0,1.0,2.0)
"

%feature("docstring") cosmology::getM4rs
"Returns the mass within 4 scale radius given the mass and concentration of a
halo

:Parameters:

-   M : Mass
-   c : concentration

:Examples:

    >>> import cosmology as cc
    >>> a = cc.cosmology(0.27, 0.73, 0.047, 0.7, 2.726, 0.82, 0.95)
    >>> a.getM4rs(M,c)
"

%include "cosmology.h"
