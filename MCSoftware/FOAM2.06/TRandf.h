#ifndef TRANDF_H
#define TRANDF_H

#include <iostream>
#include <iomanip>
#include <math.h>

#ifdef ROOT_DEF
#include "TROOT.h"
#include "TNamed.h"
#include "TObject.h"
#endif

#include"TRND.h"


/////////////////////////////////////////////////////////////////
///       RANLUX- LUXURY RANDOM NUMBER GENERATOR         ////////
/////////////////////////////////////////////////////////////////
class TRanluxEngine:  public TRND {

 private:
   int nskip, luxury;
   float float_seed_table[24];
   int i_lag,j_lag;
   float carry;
   int count24;
   const long int_modulus;
   const double mantissa_bit_24;
   const double mantissa_bit_12;
   long theSeed;                    // seed
 public:
   TRanluxEngine(void);
   TRanluxEngine(long, int);
   void SetSeedLux(long, int);
   void SetSeed(long);                  // implements RND
   double Flat(void);                   // implements RND
   void FlatArray(const int, double*);  // implements RND
   void showStatus(void);
   ///////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
   ClassDef(TRanluxEngine,1) // TRanluxEngine, ultimate r.n. of Lusher
#endif
};


//////////////////////////////////////////////////////////////////
////       RANMAR - ORDINARY RANDOM NUMBER GENERATOR      ////////
//   this will probably go...
//////////////////////////////////////////////////////////////////
class HepJamesRandom: public TRND {
 private:
  double u[97];
  double c, cd, cm;
  int i97, j97;
  long theSeed;
 public:
  HepJamesRandom();
  HepJamesRandom(long seed);
  void SetSeed(long seed);
  double Flat(void);
  void FlatArray(const int, double *);
  void ShowStatus();
  /////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(HepJamesRandom,1)  //`HepJamesRandom, RANMAR again
#endif
};


////////////////////////////////////////////////////////////////////////
////     Mersene Twistor RNG,  fast, compact, huge-period generator  ///
////////////////////////////////////////////////////////////////////////
class TRMersenneTwister : public TRND {
 private:
  unsigned int mt[624];
  int count624;
  long theSeed;
  double twoToMinus_32;
  double twoToMinus_53;
  double nearlyTwoToMinus_54;
  void powersOfTwo();
  enum{ NminusM = 227, M = 397, N = 624};
 public:
  TRMersenneTwister();
  TRMersenneTwister(long seed);
  void SetSeed(long seed);
  double Flat(void);
  void FlatArray(const int, double *);
  void ShowStatus();
  ////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(TRMersenneTwister,1) //`TRMersenneTwister, Mersene Twistor r.n. generator
#endif
};


/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//       PseudoRandom number generator TRanmarEngine (RANMAR)              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
class TRanmarEngine: public TRND {
 private:
  int       m_modcns;
  double    m_U[98], m_C,   m_CD,   m_CM,   m_twom24;
  long int  m_i97,   m_j97, m_ntot, m_ijkl, m_ntot2;
 public:
  TRanmarEngine(void);
  TRanmarEngine(long);
  TRanmarEngine(long ,long ,long );
  ~TRanmarEngine();
  void Initialize(long ,long ,long );
  void SetSeed(long);                 // implements interf. TRND
  double Flat(void);                  // implements interf. TRND
  void FlatArray(const int, double*); // implements interf. TRND
  void ShowStatus();
 private:
  //inline functions
  long int max(long int x, long int y){ if (x>=y) return(x);  else return(y);}
  long int MOD(long int x, long int y){ return(x % y);}
/////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(TRanmarEngine,1) //`TRanmarEngine, r.n. generator (RANMAR)
#endif
};

/////////////////////////////////////////////////////////////////////////////
//                End of the package TRandf                                //
/////////////////////////////////////////////////////////////////////////////
#endif

