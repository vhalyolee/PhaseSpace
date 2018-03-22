#ifndef TRND_H
#define TRND_H
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//      Class TRND                                                         //
//      Abstract class representing any type of random number generator    //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "ROOT_DEF.h"

#ifdef ROOT_DEF
#include "TROOT.h"
#include "TNamed.h"
#include "TObject.h"

class TRND: public TObject {
#else
class TRND {
#endif
 public:
  TRND() { };
  virtual ~TRND(){ };
  virtual void SetSeed(long)=0;
  virtual double Flat(void)=0;
  virtual void FlatArray(const int, double*)=0;
/////////////////////////////////////////////////////////////////////////////
#ifdef ROOT_DEF
  ClassDef(TRND,1) // Abstract class (interface) forany type of r.n. generator
#endif
};
/////////////////////////////////////////////////////////////////////////////
//                      End of     RND                                     //
/////////////////////////////////////////////////////////////////////////////
#endif
