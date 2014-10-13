#include "Physics.h"
#include "TMath.h"

class Upsilon : public Physics {
 public:
 Upsilon(TRandom* rand = 0) : Physics(rand) {
    if(_debug) cout<<"Upsilon physics running"<<endl;
    Npart = 2;
  };

  virtual void Reset(){
    for(int i = 0; i < particles.size(); ++i){
      if(particles[i]) delete particles[i];
    }
    particles.clear();

    vMax = 0.5;
    charge = 1;
    double pi = TMath::Pi();

    double M = 0.02;
    TLorentzVector upsilon;
    upsilon.SetXYZM(random->Gaus(0,vMax),random->Gaus(0,vMax),0.,M);

    double phi = random->Uniform(-pi,pi);

    TLorentzVector d1(M/2.*cos(phi),M/2.*sin(phi),0.,0.00001);
    TLorentzVector d2(-d1.X(),-d1.Y(),0.,0.00001);

    d1.Boost(upsilon.BoostVector());
    d2.Boost(upsilon.BoostVector());

    Particle* p1 = new Particle(d1.X(),d1.Y(),0.,1);
    Particle* p2 = new Particle(d2.X(),d2.Y(),0.,-1);

    Add(p1);
    Add(p2);

  }


};

