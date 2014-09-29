#ifndef Level3_h
#define Level3_h

#include "Game.h"

class Level3 : public Game {
 public:
 Level3(TCanvas* c1 = 0) :
  Game(c1)
    {
      physics->Npart = 2;
      reco->Npart = 2;
      hist = GetHist();
    }

  virtual void SetPhysics(){
    physics = new Upsilon(random);
    reco = new Reconstruction(*physics);
  }

  virtual void FillHist(){
    TLorentzVector sum;
    for(int i = 0; i < reco->recoTracks.size(); ++i){
      sum += reco->recoTracks[i]->momentum;
    }
    hist->Fill(sum.M());
    hist->Draw("");

  }

  virtual TH1D* GetHist(){
    if(hist) return hist;
    else return HistMass();
  }

};

#endif
