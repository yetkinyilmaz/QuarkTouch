#include "TVector3.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TGraph.h"

#include "TH2D.h"
#include "TH1D.h"
#include <vector>
#include <iostream>

using namespace std;

static bool _watchPropagation = 0;


class Particle{

public:
  Particle(float vx = 0.002, float vy = 0.002, float vz = 0, int c = 1) :
    position(0,0,0),
    momentum(vx,vy,vz),
    velocity(vx,vy,vz),
    charge(c)
  {;}

    void Propagate(){
    position += velocity;
    velocity += force;
  }

  TVector3 momentum;

  TVector3 position;
  TVector3 velocity;
  TVector3 force;
  int charge;

};

class Physics{
public:
  Physics(TRandom* rand = 0) : random(rand), particles(0), field(0,0,0.005) {
  };
  void Add(Particle* p){
    particles.push_back(p);
  }

  void Update(){
    for(int i = 0; i < particles.size(); ++i){
      Particle* p = particles[i];

      p->force.SetXYZ(0.001,0.002,0.0003);
      p->force = p->velocity.Cross(field)*p->charge;
      p->Propagate();

    }
  };

  vector<Particle*> particles;  
  TRandom* random;
  TVector3 field;

};

class Reconstruction : public Physics{
public:
  Reconstruction() : tracks(0) {;}
  Reconstruction(Physics& p) : 
    Physics(p), tracks(0), 
    shakeMin(0.03),shakeMax(0.04),trembleMin(0.001),trembleMax(0.01) 
  {
    for(int i = 0; i < particles.size(); ++i){

      tracks.push_back(new TGraph(0));
      recoTracks.push_back(new Particle());

      fingerAmplitudeX.push_back(random->Uniform(shakeMin, shakeMax));
      fingerFrequencyX.push_back(random->Uniform(trembleMin, trembleMax));
      fingerAmplitudeY.push_back(random->Uniform(shakeMin, shakeMax));
      fingerFrequencyY.push_back(random->Uniform(trembleMin, trembleMax));

    }
  }

  void Trace(int ip, int i, double x, double y){

    x += fingerAmplitudeX[ip]*sin(((double)i)*fingerFrequencyX[ip]);
    y += fingerAmplitudeY[ip]*sin(((double)i)*fingerFrequencyY[ip]);

    tracks[ip]->SetPoint(i,x,y);
  }

  void Fit(){
    for(int i = 0; i < tracks.size(); ++i){
      TVector3 smear(random->Gaus(0,0.01),
		     random->Gaus(0,0.01),
		     0);

      recoTracks[i]->momentum = particles[i]->momentum+smear;
    }
  }

  void Draw(TCanvas* c1 = 0, TCanvas* c2 = 0, TH1* hist = 0){

    for(int i = 0; i < tracks.size(); ++i){
      TGraph* g = new TGraph();
      g->SetLineWidth(3);
      int iii = 0;
      for(int ii = 0; ii < tracks[i]->GetN(); ii+= 30){
	double x,y;
	tracks[i]->GetPoint(ii,x,y);
	g->SetPoint(iii,x,y);
	if(c1) c1->cd();
	g->Draw("Cl");
	if(c1) c1->Update();
	++iii;
      }
      
      hist->Fill(recoTracks[i]->momentum.Mag());
      if(c2) c2->cd();
      hist->Draw("");
      if(c2) c2->Update();

    }
  };

  vector<TGraph*> tracks;
  vector<double> fingerAmplitudeX;
  vector<double> fingerFrequencyX;
  vector<double> fingerAmplitudeY;
  vector<double> fingerFrequencyY;
  vector<Particle*> recoTracks;

  double shakeMin,shakeMax,trembleMin,trembleMax;

};


void touch(){

  TH1::SetDefaultSumw2();

  TRandom* random = new TRandom();
  Physics* physics = new Physics(random);

  int Npart = 15;
  int Nstep = 500;

  double vMax = 0.01;

  int charge = 1;

  for(int i = 0; i < Npart; ++i){
    //    Particle* p = new Particle(random->Uniform(-vMax,vMax),random->Uniform(-vMax,vMax),0.,charge);
    Particle* p = new Particle(random->Gaus(0,vMax),random->Gaus(0,vMax),0.,charge);

    physics->Add(p);
    charge *= -1;
  }

  Reconstruction* reco = new Reconstruction(*physics);

  TH2D* det = new TH2D("det",";x;y",100,-1,1,100,-1,1);
  TH1D* hist = new TH1D("hist",";x;time",100,0,0.05);

  TCanvas* pad1 = new TCanvas("pad1","",800,800);
  TCanvas* pad2 = new TCanvas("pad2","",400,400);
  
  for(int i = 0; i < Nstep; ++i){
    for(int ip = 0; ip < physics->particles.size(); ++ip){
      Particle* p = (physics->particles)[ip];
      det->Fill(p->position.x(),p->position.y());
      reco->Trace(ip,i,p->position.x(),p->position.y());
    }

    if(_watchPropagation){
      pad1->cd();
      det->Draw("colz");
      pad1->Update();
    }

    physics->Update();

  }

  reco->Fit();

  pad1->cd();
  det->Draw("box");
  reco->Draw(pad1,pad2,hist);

  pad2->cd();
  hist->Draw();



}

