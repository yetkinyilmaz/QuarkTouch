#ifndef Histograms_h
#define Histograms_h

#include "TH1D.h"

TH1D* HistMass(){
  return new TH1D("hist","",100,0,0.1);
}

TH1D* HistMomentum(){
  TH1D* h = new TH1D("hist",";VELOCITY;PARTICLES",100,0,0.3);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();

  return h;
}


#endif

