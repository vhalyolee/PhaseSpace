/////////////////////////////////////////////////////////////////
//                                                             //
//       Demonstration program for  of Foam 2.x                //
//       To execute type      make Demo-run
//                            make Demo-map
//                            make Demo-run -f Makefile
//
//       Without ROOT         make DemoNR-run
//                                                             //
/////////////////////////////////////////////////////////////////
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

#include <math.h>

#include "TRND.h"
#include "TRandf.h"
#include "TFDISTR.h"
#include "TFOAM_INTEGRAND.h"
#include "TFOAM.h"        // ROOT headers hidden here
 
#ifdef ROOT_DEF
#include "TFile.h"
#else
#include "TFHST.h"
#endif


#define SP12 setprecision(7) << setw(12)
#define SW18 setw(18)

#ifdef ROOT_DEF
TFile RootFile("rdemo.root","RECREATE","histograms");
#endif

int main(void){
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//        MAIN PROGRAM starts here                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
  long   loop;
  double MCresult,MCerror,MCwt;
  int    FunType;
//=========================================================
  long NevTot   = 1000000;   // Total MC statistics
  int  nDim     =       2;   // SIMPLICAL   subspace
  int  kDim     =       0;   // HYP-CUBICAL subspace
  int  TotDim   = kDim+ nDim;// total dimension
  int  nCells   =    5000;   // Number of Cells
  int  nSampl   =   20000;   // Number of MC events per cell in build-up
  int  nBin     =       8;   // Number of bins in build-up
  int  OptRej   =       1;   // =0, weighted events;  =1, wt=1 events
  int  OptDrive =       2;   // (D=2) Option, type of Drive =0,1,2 for TrueVol,Sigma,WtMax
  int  OptEdge  =       0;   // (D=0) Option, vertices are included in the sampling (1) or not (0)
  int  OptOrd   =       0;   // (D=0) Option, single simplex (1) or nDim! simplices (0)
  int  OptPeek  =       0;   // (D=0) Option, choice of cell in buid-up maximum(0), random(1)
  int  OptMCell =       1;   // (D=1) Option, Mega-Cell = slim memory
  int  OptVert  =       1;   // (D=1) Vertices are not stored
  int  EvPerBin =      25;   // Maximum events (equiv.) per bin in buid-up
  double MaxWtRej =   1.0;   // Maximum wt for rejection, for OptRej=1
  int  Chat     =       1;   // Chat level
//=========================================================
  FunType =10; // see TFDISTR.cxx
  TFDISTR *Density1 = new TFDISTR(FunType); // Create integrand function
  TRND *PseRan   = new TRanmarEngine();  // Create random number generator
  TFOAM   *FoamX    = new TFOAM("FoamX");   // Create Simulator
//=========================================================
  double pifac = 4.0*atan(1.0);
  cout<<"*****   Demonstration Program for Foam version "<<FoamX->GetVersion()<<"    *****"<<endl;
  FoamX->SetnDim(        nDim);
  FoamX->SetkDim(        kDim);
  FoamX->SetnCells(      nCells);
  FoamX->SetnSampl(      nSampl);
  FoamX->SetnBin(        nBin);
  FoamX->SetOptRej(      OptRej);
  FoamX->SetOptDrive(    OptDrive);
  FoamX->SetOptPeek(     OptPeek);
  FoamX->SetOptEdge(     OptEdge);
  FoamX->SetOptOrd(      OptOrd);
  FoamX->SetOptMCell(    OptMCell);
  FoamX->SetOptVert(     OptVert);
  FoamX->SetEvPerBin(    EvPerBin);
  FoamX->SetMaxWtRej(    MaxWtRej);
  FoamX->SetChat(        Chat);
//===============================
  FoamX->Initialize(PseRan,  Density1 ); // Initialize simulator
#ifdef ROOT_DEF
  FoamX->Write("FoamX");    // Writing Foam on the disk, TESTING PERSISTENCY!!!
#endif
//===============================
  long nCalls=FoamX->GetnCalls();
  if((nDim==2)||(kDim==2)){ FoamX->LaTexPlot2dim("Demo-map.tex"); }
  cout << "====== Initialization done, entering MC loop" << endl;
  double XCrude = FoamX->GetPrimary();
//======================================================================
//======================================================================
  //cout<<" About to start MC loop: ";  cin.getline(question,20);
  double *MCvect =new double[TotDim]; // vector generated in the MC run
//======================================================================
#ifdef ROOT_DEF
  TH1D  *hst_Wt = new TH1D("hst_Wt" ,  "Main weight of Foam",25,0,1.25);
  hst_Wt->Sumw2();
#else
  TFHST *hst_Wt = new TFHST(0.0,1.25, 25);
#endif
//======================================================================
  for(loop=0; loop<NevTot; loop++){
//===============================
    FoamX->MakeEvent();           // generate MC event
//===============================
    FoamX->GetMCvect( MCvect);
    MCwt=FoamX->GetMCwt();
    hst_Wt->Fill(MCwt,1.0);    
    double Xvec[nDim];
    FoamX->GetMCvect(Xvec);
    double myWt = 1.0;
    double myX = Xvec[0];
    if ((myX>0.0) && (myX<1.0)) {
          double s1 = sin(pifac*myX);
          double s2 = 4.0*myX*(1.0-myX);
          myWt = s1*s1/s2/s2;
    }
    MCwt *= myWt;
    cout<<"XWt "<<myX<<"  "<<MCwt<<endl;
    if(loop<15){
      cout<<"MCwt= "<<SP12<<MCwt<<",  ";
      cout<<"MCvect= ";
      for ( int k=0 ; k<TotDim ; k++) cout<<SP12<<MCvect[k]<<" "; cout<< endl;
    }
    if( ((loop)%100000)==0 ){
      cout<<"   loop= "<<loop<<endl;
    }
  }
//======================================================================
//======================================================================
  cout << "====== Events generated, entering Finalize" << endl;

  hst_Wt->Print("all");
  double eps = 0.0005;
  double Effic, WtMax, AveWt, Sigma;
  double IntNorm, Errel;
  FoamX->Finalize(   IntNorm, Errel);     // final printout
  FoamX->GetIntegMC( MCresult, MCerror);  // get MC integral
  FoamX->GetWtParams(eps, AveWt, WtMax, Sigma); // get MC wt parameters
  Effic=0; if(WtMax>0) Effic=AveWt/WtMax;
  cout << "================================================================" << endl;
  cout << " MCresult= " << MCresult << " +- " << MCerror << " RelErr= "<< MCerror/MCresult << endl;
  cout << " Dispersion/<wt>= " << Sigma << endl;
  cout << "      <wt>/WtMax= " << Effic <<",    for epsilon = "<<eps << endl;
  cout << " nCalls (initialization only) =   " << nCalls << endl;
  cout << " XCrude         = " << XCrude << endl;
  cout << "================================================================" << endl;
// Analytic results for Herwig case:
    double c1 = 3.2168353;
    double svar = 7000.0;
    svar *= svar;
    double z1 = 60.0/7000.0;
    z1 *= z1;
    double z2 = 110.0/7000.0;
    z2 *= z2;
// Dist 8 or 9
//    double integral = (c1/svar)*((log(z2)+1.0)/z2 - (log(z1)+1.0)/z1);
// Dist 10
      double integral = 5.0/8.0 - cos(2.0)/4.0 + cos(4.0)/8.0;
    cout << "Analytic result = " << integral <<endl;

  delete [] MCvect;
  //
#ifdef ROOT_DEF
  RootFile.ls();
  RootFile.Write();
  RootFile.Close();
#endif
  cout << "***** End of Demonstration Program  *****" << endl;
} // end of Main


