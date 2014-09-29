#ifndef Reconstruction_h
#define Reconstruction_h

#include "Physics.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1D.h"
#include <iostream>

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
      TLorentzVector smear(random->Gaus(0,0.0000000000001),
			   random->Gaus(0,0.0000000000001),
			   0,
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

      if(hist){
        hist->Fill(recoTracks[i]->momentum.P());
        if(c2) c2->cd();
        hist->Draw("");
        if(c2) c2->Update();
      }

    }

  };
  virtual void Reset(std::vector<Particle*> p){

    particles = p;
    for(int i = 0; i < tracks.size(); ++i){
      if(_debug)std::cout<<"deleting track "<<i<<std::endl;
      if(tracks[i]) delete tracks[i];
      if(_debug)std::cout<<"deleting recotrack "<<i<<std::endl;
      if(recoTracks[i]) delete recoTracks[i];
    }

    if(_debug)std::cout<<"pointers deleted"<<std::endl;

    tracks.clear();
    recoTracks.clear();

    fingerAmplitudeX.clear();
    fingerFrequencyX.clear();
    fingerAmplitudeY.clear();
    fingerFrequencyY.clear();


    if(_debug)std::cout<<"Creating reco objects"<<std::endl;
    for(int i = 0; i < particles.size(); ++i){
      if(_debug)std::cout<<"object "<<i<<std::endl;
      tracks.push_back(new TGraph(0));
      recoTracks.push_back(new Particle());

      fingerAmplitudeX.push_back(random->Uniform(shakeMin, shakeMax));
      fingerFrequencyX.push_back(random->Uniform(trembleMin, trembleMax));
      fingerAmplitudeY.push_back(random->Uniform(shakeMin, shakeMax));
      fingerFrequencyY.push_back(random->Uniform(trembleMin, trembleMax));
    }
  }

  std::vector<TGraph*> tracks;
  std::vector<double> fingerAmplitudeX;
  std::vector<double> fingerFrequencyX;
  std::vector<double> fingerAmplitudeY;
  std::vector<double> fingerFrequencyY;
  std::vector<Particle*> recoTracks;

  double shakeMin,shakeMax,trembleMin,trembleMax;

};

#endif
