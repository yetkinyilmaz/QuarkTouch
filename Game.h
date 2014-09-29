#ifndef Game_h
#define Game_h

#include "Globals.h"
#include "Physics.h"
#include "Reconstruction.h"
#include "Histograms.h"

#include "TH2D.h"

class Game {
 public:

 Game(TCanvas* c1 = 0) : hist(0) {
    det = new TH2D("det",";x;y",100,-1,1,100,-1,1);
    pad1 = c1;
    Nstep = 200;
    random = new TRandom();
    SetPhysics();
  }

  virtual void SetPhysics(){
    physics = new Physics(random);
    reco = new Reconstruction(*physics);
  }

  void Generate(){

    det->Reset();
    physics->Reset();
    reco->Reset(physics->particles);

    if(_debug)cout<<"Particles created"<<endl;

    for(int i = 0; i < Nstep; ++i){
      for(int ip = 0; ip < physics->particles.size(); ++ip){
        if(_debug)cout<<"Propagating particle "<<ip<<endl;
        Particle* p = (physics->particles)[ip];
        det->Fill(p->position.x(),p->position.y());
        reco->Trace(ip,i,p->position.x(),p->position.y());
      }

      if(_watchPropagation){
        pad1->cd();
        det->Draw("box");
        pad1->Update();
      }

      physics->Update();

    }
    if(_debug)cout<<"Particles propagated"<<endl;

    pad1->cd();
    det->Draw("box");
    reco->Fit();
  }

  virtual void FillHist(){
    return;
  }

  virtual TH1D* GetHist(){return 0;}

  TRandom* random;
  Physics* physics;
  Reconstruction* reco;
  TH2D* det;
  TH1D* hist;
  TCanvas* pad1;
  int Nstep;

};

#endif
