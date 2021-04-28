
#ifndef B1ActionInitialization_h
#define B1ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "B1DetectorConstruction.hh"
#include "ParamsSingleton.hh"

class B1ActionInitialization : public G4VUserActionInitialization
{
  public:
    B1ActionInitialization(B1DetectorConstruction* geometry);
    virtual ~B1ActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    B1DetectorConstruction* fdet;
   	ParamsSingleton* pp = ParamsSingleton::Instance();

};


#endif

    
