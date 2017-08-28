// A very simple script to make a plot showing the distribution of random
// numbers in a 2-D plane using ROOT's TRandom3 (Mersenne twister) generator.

// This script is designed to be compiled using the included Makefile.
// But you can also run it interactively within ROOT (using .x) if you change
// the name of main to "PlotMersenneTwister".

#include "TRandom3.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TStyle.h"
#include "TMarker.h"

const int npoints = 1000;

int main(void) {
  // set up PRNG. 0 will automatically generate a new seed using TUUID
  // or if you want a reproducible plot, use a fixed seed
  TRandom3 *r = new TRandom3(0);

  // set up canvas and dummy histogram
  TCanvas *c = new TCanvas("c", "c", 600, 600);
  TH2F *h = new TH2F("h", "h", 5, 0, 1, 5, 0, 1);

  gStyle->SetOptStat(0);
  h->Draw();
  h->SetTitle("2-D Mersenne Twister sequence");
  h->GetXaxis()->SetTitle("x");
  h->GetYaxis()->SetTitle("y");

  // Generate points and plot them on the histogram
  double point[2];
  for (int i=0; i<npoints; ++i) {
    r->RndmArray(2, point);

    TMarker *m = new TMarker(point[0], point[1], kPlus);
    m->Draw();
  }

  c->Print("MersenneTwister.png");
}
