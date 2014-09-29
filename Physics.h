#ifndef Physics_h
#define Physics_h


#include "Particle.h"
#include "TRandom.h"
#include <vector>

class Physics{
 public:
 Physics(TRandom* rand = 0) : random(rand), particles(0), field(0,0,0.005) {
    Npart = 10;
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

    vMax = 0.01;
    charge = 1;
    for(int i = 0; i < Npart; ++i){
      Particle* p = ProduceParticle();
      Add(p);
      charge *= -1;
    }

  }

  virtual Particle* ProduceParticle(){
    return new Particle(random->Gaus(0,vMax),random->Gaus(0,vMax),0.,charge);
  }

  std::vector<Particle*> particles;
  TRandom* random;
  TVector3 field;
  int Npart;
  int charge;
  double vMax;
};

#endif
