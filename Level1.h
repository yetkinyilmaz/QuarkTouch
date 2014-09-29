#ifndef Level1_h
#define Level1_h

#include "Game.h"

class Level1 : public Game {
 public:
 Level1(TCanvas* c1 = 0) :
  Game(c1)
    {
      physics->Npart = 1;
      reco->Npart = 1;
      hist = GetHist();
    }

  virtual TH1D* GetHist(){
    std::cout<<"Momentum histogram"<<endl;
    if(hist) return hist;
    else return HistMomentum();
  }

  virtual void FillHist(){

    double ze = 0;
    for(int i = 0; i < reco->recoTracks.size(); ++i){
      ze += reco->recoTracks[i]->momentum.P();
      hist->Fill(ze);
    }
    hist->Draw("");
  }

};

#endif
