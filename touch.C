#include "TVector3.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TGraph.h"

#include "TH2D.h"
#include "TH1D.h"
#include <vector>
#include <iostream>

using namespace std;

static bool _watchPropagation = 1;
static bool _debug = 0;
static int _drawSpeed = 20;

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

  virtual void Reset(){
    for(int i = 0; i < particles.size(); ++i){
      if(particles[i]) delete particles[i];
    }
    particles.clear();

    int Npart = 10;

    double vMax = 0.01;
    int charge = 1;
    for(int i = 0; i < Npart; ++i){
      Particle* p = new Particle(random->Gaus(0,vMax),random->Gaus(0,vMax),0.,charge);
      Add(p);
      charge *= -1;
    }


  }

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
  {;}

  void Trace(int ip, int i, double x, double y){
    x += fingerAmplitudeX[ip]*sin(((double)i)*fingerFrequencyX[ip]);
    y += fingerAmplitudeY[ip]*sin(((double)i)*fingerFrequencyY[ip]);
    tracks[ip]->SetPoint(i,x,y);
  }

  void Fit(){
    for(int i = 0; i < tracks.size(); ++i){
      TVector3 smear(random->Gaus(0,0.00000001),
		     random->Gaus(0,0.00000001),
		     0);

      recoTracks[i]->momentum = particles[i]->momentum+smear;
    }
  }

  void Draw(TCanvas* c1 = 0, TCanvas* c2 = 0, TH1* hist = 0){

    for(int i = 0; i < tracks.size(); ++i){
      TGraph* g = new TGraph();
      g->SetLineWidth(3);
      int iii = 0;
      for(int ii = 0; ii < tracks[i]->GetN(); ii+= _drawSpeed){
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

  virtual void Reset(vector<Particle*> p){

    particles = p;
    for(int i = 0; i < tracks.size(); ++i){
      if(_debug)cout<<"deleting track "<<i<<endl;
      if(tracks[i]) delete tracks[i];
      if(_debug)cout<<"deleting recotrack "<<i<<endl;
      if(recoTracks[i]) delete recoTracks[i];
    }

    if(_debug)cout<<"pointers deleted"<<endl;    

    tracks.clear();
    recoTracks.clear();

    fingerAmplitudeX.clear();
    fingerFrequencyX.clear();
    fingerAmplitudeY.clear();
    fingerFrequencyY.clear();


    if(_debug)cout<<"Creating reco objects"<<endl;
    for(int i = 0; i < particles.size(); ++i){
      if(_debug)cout<<"object "<<i<<endl;
      tracks.push_back(new TGraph(0));
      recoTracks.push_back(new Particle());

      fingerAmplitudeX.push_back(random->Uniform(shakeMin, shakeMax));
      fingerFrequencyX.push_back(random->Uniform(trembleMin, trembleMax));
      fingerAmplitudeY.push_back(random->Uniform(shakeMin, shakeMax));
      fingerFrequencyY.push_back(random->Uniform(trembleMin, trembleMax));
    }
  }


  vector<TGraph*> tracks;
  vector<double> fingerAmplitudeX;
  vector<double> fingerFrequencyX;
  vector<double> fingerAmplitudeY;
  vector<double> fingerFrequencyY;
  vector<Particle*> recoTracks;

  double shakeMin,shakeMax,trembleMin,trembleMax;

};


class Game{
public:

  Game(TH1* h = 0, TCanvas* c1 = 0){
    random = new TRandom();
    physics = new Physics(random);
    reco = new Reconstruction(*physics);
    det = new TH2D("det",";x;y",100,-1,1,100,-1,1);
    hist = h;
    pad1 = c1;
  }

  void Generate(){

    det->Reset();
    physics->Reset();
    reco->Reset(physics->particles);

    int Nstep = 300;

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

  TRandom* random;
  Physics* physics;
  Reconstruction* reco;
  TH2D* det;
  TH1* hist;
  TCanvas* pad1;
};

void touch(){

  TH1::SetDefaultSumw2();

  int Nevents = 200;

  TCanvas* pad1 = new TCanvas("pad1","",800,800);
  TCanvas* pad2 = new TCanvas("pad2","",400,400);
  TH1D* hist = new TH1D("hist",";x;time",100,0,0.05);

  Game* game = new Game(hist,pad1);

  for(int i = 0; i< Nevents; ++i){
    if(_debug)cout<<"Event : "<<i<<endl;
    game->Generate();
    game->reco->Draw(pad1,pad2,hist);
  }


}



