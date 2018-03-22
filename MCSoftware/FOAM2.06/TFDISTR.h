#ifndef TFDISTR_H
#define TFDISTR_H
#include "ROOT_DEF.h"

#ifdef ROOT_DEF
#include "TROOT.h"
#endif
#include "TFOAM_INTEGRAND.h"

class TFOAM_INTEGRAND;
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//      Class of testing functions                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
class TFDISTR: public TFOAM_INTEGRAND {
  private:
    int    m_Type;                      // Type of function
    double m_Pi;                        // constant pi
  public:
    TFDISTR(int t=0);                     // Constructor
    virtual ~TFDISTR(void);               // Destructor
// methods
  public:
    double Density(int, double *);      // MainTesting functions
    double Sphere( int, double *);
    double Void(   int, double *);
    double Diagon( int, double *);
    double Camel(  int, double *);
    double Gauss(  int, double *);
    double Ridge(  int, double *);
    double Test1d(  int, double *);
    double Herwig1( int, double *);
    double Herwig2( int, double *);
////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
    ClassDef(TFDISTR,1) //Class of testing functions for FOAM
#endif
};
#endif
