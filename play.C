

#include "Level1.h"
#include "Level2.h"
#include "Level3.h"

#include "Histograms.h"


void play(){

  TH1::SetDefaultSumw2();

  int Nevents = 5;
  int level = 1;

  TCanvas* pad1 = new TCanvas("pad1","",800,800);
  TCanvas* pad2 = new TCanvas("pad2","",400,400);
  TH1D* hist = 0;

  Level2* game = new Level2(pad1);

  /*
  switch(level) {
  case 1 : game = new Level1(pad1);
  //    case 2 : game = new Level2(pad1);
  //    case 3 : game = new Level3(pad1);
  }
  */

  game->SetPhysics();
  hist = game->GetHist();

  for(int i = 0; i< Nevents; ++i){
    if(_debug)cout<<"Event : "<<i<<endl;
    game->Generate();
    if(0){
      game->reco->Draw(pad1,pad2,hist);
    }else{
      if(!_skipToResult) game->reco->Draw(pad1,pad2,0);
      pad2->cd();
      game->FillHist();
      if(!_skipToResult) pad2->Update();
    }
  }


}



