#include "Physics.h"

class Upsilon : public Physics {
 public:
 Upsilon(TRandom* rand = 0) : Physics(rand) {
    if(_debug) cout<<"Upsilon physics running"<<endl;
    Npart = 2;
  };

};

