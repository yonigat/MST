
#ifndef B1BOptrFDTrackData_hh
#define B1BOptrFDTrackData_hh

class B1BOptrFD;
#include "G4VAuxiliaryTrackInformation.hh"

enum class FDState { normal, toBeFreeFlight };

class B1BOptrFDTrackData : public G4VAuxiliaryTrackInformation {

friend class B1BOptrFD;

public:
  B1BOptrFDTrackData(  const B1BOptrFD* optr );
  ~B1BOptrFDTrackData();

  // -- from base class:
  void Print() const;

  // -- Get methods:
  G4bool  IsFreeFromBiasing() const
  { return ( fFDState == FDState::normal);}
  G4bool IsPrimary()
  {return primary;}
  void SetSecondary()
  {primary=false;}
  const B1BOptrFD* GetOptr()
  {return fFDOperator;}

  FDState             fFDState;


  // -- no set methods are provided : sets are made under exclusive control of B1BOptrFS objects through friendness.

private:
  const B1BOptrFD*    fFDOperator;
  G4bool              primary;

  void Reset()
    {
      fFDOperator = nullptr;
      fFDState    = FDState::normal;
    }



};

#endif
