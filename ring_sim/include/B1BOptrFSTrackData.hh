
#ifndef B1BOptrFSTrackData_hh
#define B1BOptrFSTrackData_hh

class B1BOptrFS;
#include "G4VAuxiliaryTrackInformation.hh"

enum class FSState { start, finish, toBeFreeFlight, toBeSplitCompt, toBeSplitRayl };

class B1BOptrFSTrackData : public G4VAuxiliaryTrackInformation {

friend class B1BOptrFS;

public:
  B1BOptrFSTrackData(  const B1BOptrFS* optr );
  ~B1BOptrFSTrackData();

  // -- from base class:
  void Print() const;

  // -- Get methods:
  G4bool  IsFreeFromBiasing() const
  { return ( fFSState == FSState::start);}
  G4bool IsPrimary()
  {return primary;}
  void SetSecondary()
  {primary=false;}
  const B1BOptrFS* GetOptr()
  {return fFSOperator;}

  FSState             fFSState;


  // -- no set methods are provided : sets are made under exclusive control of B1BOptrFS objects through friendness.

private:
  const B1BOptrFS*    fFSOperator;
  G4bool              primary;

  void Reset()
    {
      fFSOperator = nullptr;
      fFSState    = FSState::start;
    }



};

#endif
