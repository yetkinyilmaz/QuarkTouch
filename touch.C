#include "TVector3.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include <vector>
#include <iostream>

using namespace std;

class Particle{


public:
  Particle() :
    position(-0.5,0.5,0),
    velocity(0.002,0.002,0)
  {;}

    void Propagate(){
    position += velocity;
    velocity += force;
  }


  TVector3 position;
  TVector3 velocity;
  TVector3 force;


};

class Physics{
public:
  Physics() : particles(0), field(0,0,0.005) {
  };
  void Add(Particle* p){
    particles.push_back(p);
  }

  void Update(){
    for(int i = 0; i < particles.size(); ++i){
      //      particles[i]->force.SetXYZ(0.001,0.002,0.0003);
      particles[i]->force = particles[i]->velocity.Cross(field);
      particles[i]->Propagate();

    }
  };

  vector<Particle*> particles;  
  TVector3 field;

};


void touch(){

  TH1::SetDefaultSumw2();

  Particle* p = new Particle();
  Physics* physics = new Physics();

  physics->Add(p);

  TH2D* det = new TH2D("det",";x;y",50,-1,1,50,-1,1);
  TH1D* hist = new TH1D("hist",";x;time",100,-2,2);

  int N = 10000;

  TCanvas* pad1 = new TCanvas("pad1","",800,800);
  TCanvas* pad2 = new TCanvas("pad2","",400,400);
  
  for(int i = 0; i < N; ++i){
    physics->Update();
    det->Fill(p->position.x(),p->position.y());
    hist->Fill(p->position.x());


    pad1->cd();
    det->Draw("colz");
    pad1->Update();

    pad2->cd();
    hist->Draw();
    pad2->Update();
    physics->Update();

  }

  pad2->cd();
  hist->Draw();

}

