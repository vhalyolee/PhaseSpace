#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <math.h>
#include "TFDISTR.h"

///////////////////////////////////////////////////////////////////////////////
extern "C" double circe_ee_(const double&, const double&);

#ifdef ROOT_DEF
ClassImp(TFDISTR)
#endif
///////////////////////////////////////////////////////////////////////////////
//                          Class TFDISTR                                    //
//                                                                           //
// Class TFDISTR is collection of the testing functions                      //
// Inherits  from abstract class TFOAM_INTEGRAND!                            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
TFDISTR::TFDISTR(int typ){
  m_Type=typ;      // type of integrand
  m_Pi = 3.1415926535897932e0;
}
//-----------------------------------------------------------------------------
TFDISTR::~TFDISTR(void){ }

//-----------------------------------------------------------------------------
double TFDISTR::Density(int nDim, double *Xarg){
  double rho;
  rho = 1e-100;
  if (      m_Type == 1) 
    rho = Sphere( nDim,Xarg);
  else if ( m_Type == 2) 
    rho = Void(   nDim,Xarg);
  else if ( m_Type == 3) 
    rho = Diagon( nDim,Xarg);
  else if ( m_Type == 4) 
    rho = Camel( nDim,Xarg);
  else if ( m_Type == 5) 
    rho = Gauss( nDim,Xarg);
  else if ( m_Type == 6) 
    rho = circe_ee_(Xarg[0], Xarg[1]);
  else if ( m_Type == 7) 
    rho = Ridge( nDim,Xarg);
  else if ( m_Type == 8) 
    rho = Herwig1( nDim,Xarg);
  else if ( m_Type == 9) 
    rho = Herwig2( nDim,Xarg);
  else if ( m_Type == 10) 
    rho = Test1d( nDim,Xarg);
  else{
    cout<<" TFDISTR::Density: !!!!!! WRONG m_Type= "<<m_Type<<endl;
    exit(0);
  }
  return rho;
}
//-----------------------------------------------------------------------------
// Testing function of Peter Lepage
//-----------------------------------------------------------------------------
double TFDISTR::Camel(int nDim, double *Xarg){
  double Fun1,Fun2,R1,R2;
  double pos1=1e0/3e0;
  double pos2=2e0/3e0;
  double Gam1= 0.100e0;  // as in JPC
  double Gam2= 0.100e0;  // as in JPC
  double sPi = sqrt(m_Pi);
  double xn1=1e0;
  double xn2=1e0;
  int i;
  R1=0;
  R2=0;
  for(i = 0 ; i<nDim ; i++){
    R1=R1+(Xarg[i] -pos1)*(Xarg[i] -pos1);
    R2=R2+(Xarg[i] -pos2)*(Xarg[i] -pos2);
    xn1=xn1*Gam1*sPi;
    xn2=xn2*Gam2*sPi;      
  }
  R1   = sqrt(R1);
  R2   = sqrt(R2);
  Fun1 = exp(-(R1*R1)/(Gam1*Gam1))/xn1;  // Gaussian delta-like profile
  Fun2 = exp(-(R2*R2)/(Gam2*Gam2))/xn2;  // Gaussian delta-like profile
  return 0.5e0*(Fun1+ Fun2);
}
//-----------------------------------------------------------------------------
// Testing function Gauss
//-----------------------------------------------------------------------------
double TFDISTR::Gauss(int nDim, double *Xarg){
  double Fun1,R1;
  double pos1=1e0/3e0;
  double Gam1= 0.100e0;  // as in Camel
  double sPi = sqrt(m_Pi);
  double xn1=1e0;
  int i;
  R1=0;
  for(i = 0 ; i<nDim ; i++){
    R1=R1+(Xarg[i] -pos1)*(Xarg[i] -pos1);
    xn1=xn1*Gam1*sPi;
  }
  R1   = sqrt(R1);
  Fun1 = exp(-(R1*R1)/(Gam1*Gam1))/xn1;  // Gaussian delta-like profile
  return Fun1;
}
//-----------------------------------------------------------------------------
//  n-dimensional thin sphere
//-----------------------------------------------------------------------------
double TFDISTR::Sphere(int nDim, double *Xarg){
  double f1,f2;
  double *Pos = new double[nDim];
  double Radius=0.35, Gam=0.02;
  double R,r;
  int i;
  for(i = 0 ; i<nDim ; i++){
    if( (2*(i/2)-i) == 0 )  
      Pos[i]=0.25;   // even
    else
      Pos[i]=0.40;}  // odd
  R=0;
  for(i = 0 ; i<nDim ; i++){
    r=Xarg[i]-Pos[i];
    R= R+r*r;}
  R=sqrt(R);
  f1 = Gam/((R-Radius)*(R-Radius) + Gam*Gam)/sqrt(m_Pi); // B-W profile
  f2 = exp(-R*R/(Gam*Gam))/Gam/sqrt(m_Pi);               // Gaussian profile
  delete [] Pos;
  return f1;
}
//-----------------------------------------------------------------------------
//  diagonal in n-dimensions
//-----------------------------------------------------------------------------
double TFDISTR::Diagon(int nDim, double *Xarg){
  double f1, Gam = 0.05e0;
  double *vtr = new double[nDim];
  double sum,R,r;
  int i;
  sum=0;
  for(i = 0 ; i<nDim ; i++)
    sum=sum+Xarg[i];
  R=0;
  for(i = 0 ; i<nDim ; i++){
    r=Xarg[i]-sum/nDim;
    R= R+r*r;}
  R=sqrt(R);
  f1 = Gam/(R*R + Gam*Gam)/sqrt(m_Pi); // B-W profile
  delete [] vtr;
  return f1;
}
//-----------------------------------------------------------------------------
//  simple 1d test: integral over unit interval = 1/2 
//----------------------------------------------------------------------------
double TFDISTR::Test1d(int nDim, double *Xarg){
//   double s1 = sin(m_Pi*Xarg[0]);
//    double f1 = s1*s1;
    double f1 = 4.0*Xarg[0]*(1.0-Xarg[0]);
    return f1*f1;
}
//-----------------------------------------------------------------------------
//  surface of n-dimensional hypercubic
//-----------------------------------------------------------------------------
double TFDISTR::Void(int nDim, double *Xarg){
  double f1;
  double  Gam = 0.05;
  int i;
  f1 = 0e0;
  for(i = 0 ; i<nDim ; i++){
    if(      Xarg[i] <Gam )  f1 = 1/Gam;
    if( (1.0-Xarg[i])<Gam )  f1 = 1/Gam;
  }
  return f1;
}
//-----------------------------------------------------------------------------
//  2-dim sharp ridge of Perret-Gallix
//-----------------------------------------------------------------------------
double TFDISTR::Ridge(int nDim, double *Xarg){
  double f1;
  double  alf;
  alf = 0.01;
  alf = 0.0001;
  alf = 1e-6;
  double D= Xarg[0]+Xarg[1]-1;
  f1 = 2*alf*Xarg[1]/( D*D  +alf*alf);
  return f1;
}
//----------------------------------------------------------------------------
// Herwig style integrand : generates z, x1
double TFDISTR::Herwig1(int nDim, double *Xarg) {
    double pi = 3.1415926535897932;
    double gnanob = 389379.66;
    double alfinv = 137.0359895;
    double c1 = 3.2168353;
    double svar = 7000.0;
    svar *= svar;
    double sig0nb = 4.0*pi*gnanob/(alfinv*alfinv*3.0*svar);
    double z1 = 60.0/7000.0;
    z1 *= z1;
    double z2 = 110.0/7000.0;
    z2 *= z2;
    double r1 = Xarg[0];
//  Version 1
    double z = z1*(1.0 - r1) + z2*r1;
    double jac1 = z2 - z1;
//  Version 2
//  z = 1.0/((1.0-r1)/z1 + r1/z2);
//  jac1 = z*z*(1.0/z1 - 1.0/z2);

    double r2 = Xarg[1];
    double x1 = pow(z, r2);
    double jac2 = -log(z); 
    double f1 = jac1*jac2*c1/svar/z/z;
//    double f1 = jac1*jac2*sig0nb/z;
    return f1;
}
// Herwig style integrand : generates x1, x2
double TFDISTR::Herwig2(int nDim, double *Xarg) {
    double f1 = 0.;
    double c1 = 3.2168353;
    double svar = 7000.0;
    svar *= svar;
    double z1 = 60.0/7000.0;
    z1 *= z1;
    double z2 = 110.0/7000.0;
    z2 *= z2;
    double r1 = Xarg[0];
    double x1 = r1 + (1.0-r1)*z1;
    double x2max = 1.0;
    double x2min = z1/x1;
    if (x1 > z2) {
        x2max = z2/x1;
     }
    double r2 = Xarg[1];
    double x2 = (1.0 - r2)*x2min + r2*x2max;
    double z = x1*x2;
    double jac = (1.0 - z1)*(x2max - x2min);
    if ((z > z1) && (z < z2)) {
       f1 = jac*c1/(svar*z)/z;
    } else {
       cout <<"Error: " <<  x1 << " " << x2 << " " << z << endl;
    }
    return f1;
}
///////////////////////////////////////////////////////////////////////////////
//           end of Class TFDISTR                                            //
///////////////////////////////////////////////////////////////////////////////


