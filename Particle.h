#ifndef Particle_h
#define Particle_h

#include "TLorentzVector.h"
#include "TVector3.h"

class Particle{

 public:
 Particle(float vx = 0.002, float vy = 0.002, float vz = 0, int c = 1) :
  position(0,0,0),
    velocity(vx,vy,vz),
    charge(c)
    {
      momentum = TLorentzVector(vx,vy,vz,velocity.Mag());
    }

  void Propagate(){
    position += velocity;
    velocity += force;
  }

  TLorentzVector momentum;

  TVector3 position;
  TVector3 velocity;
  TVector3 force;
  int charge;

};

#endif
