

#include "Level2.h"
#include "Histograms.h"


void touch(){

  TH1::SetDefaultSumw2();

  int Nevents = 3000;

  TCanvas* pad1 = new TCanvas("pad1","",800,800);
  TCanvas* pad2 = new TCanvas("pad2","",400,400);
  TH1D* hist = HistMass();;

  Game* game = new Level2(hist,pad1);
  game->SetPhysics();

  for(int i = 0; i< Nevents; ++i){
    if(_debug)cout<<"Event : "<<i<<endl;
    game->Generate();
    if(0){
      game->reco->Draw(pad1,pad2,hist);
    }else{
      if(!_skipToResult) game->reco->Draw(pad1,pad2,0);
      pad2->cd();
      game->FillHist(hist);
      if(!_skipToResult) pad2->Update();
    }
  }


}



