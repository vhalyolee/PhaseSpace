// Testing persistency
// make TestPers
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
using namespace std; 

#ifdef ROOT_DEF
#include "TROOT.h"
#include "TFile.h"
#endif

#include "TRND.h"
#include "TFOAM_INTEGRAND.h"
#include "TFOAM.h"
#include "TFDISTR.h"

#define SP12 setprecision(7) << setw(12)

main(int argc, char **argv)
{
  cout<<"====================== TestVector ================================"<<endl;
  TFile fileA("rdemo.root");
  fileA.cd();
  cout<<"------------------------------------------------------------------"<<endl;
  fileA.ls();
  cout<<"------------------------------------------------------------------"<<endl;
  fileA.Map();
  //cout<<"------------------------------------------------------------------"<<endl;
  //fileA.ShowStreamerInfo();
  cout<<"------------------------------------------------------------------"<<endl;
  fileA.GetListOfKeys()->Print();
  cout<<"------------------------------------------------------------------"<<endl;
  //*******************************************
  TFOAM  *FoamX = (TFOAM*)fileA.Get("FoamX");
  FoamX->LinkCells();
  //*******************************************
  //FoamX->PrintCells();
  //FoamX->PrintVertices();
  FoamX->CheckAll(1);
  cout<<"------------------------- MC loop will start----------------------"<<endl;
  int TotDim = FoamX->GetTotDim();
  double *MCvect =new double[ TotDim ]; // vector generated in the MC run
  int Nevt =15;
  for(int loop=0; loop< Nevt; loop++){
    FoamX->MakeEvent();                 // generate MC event
    FoamX->GetMCvect( MCvect);
    double MCwt=FoamX->GetMCwt();
    cout<<"MCwt= "<<SP12<<MCwt<<",  ";
    cout<<"MCvect= ";
    for ( int k=0 ; k<TotDim ; k++) cout<<SP12<<MCvect[k]<<" "; cout<< endl;
  }
  TRND  *RNGen =   FoamX->GetPseRan();        // get pointer of restored RN generator
  TFDISTR  *RHO   =   (TFDISTR*)FoamX->GetRho(); // get pointer of restored distribution
  //
  cout<<"===================== TestPers FINISHED ======================="<<endl;
  return 0;
}
//_____________________________________________________________________________
