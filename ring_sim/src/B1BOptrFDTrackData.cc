//
#include "B1BOptrFDTrackData.hh"
#include "B1BOptrFD.hh"


B1BOptrFDTrackData::B1BOptrFDTrackData( const B1BOptrFD* optr )
: G4VAuxiliaryTrackInformation(),
  fFDOperator( optr ),
  primary(true)
{
  fFDState = FDState::normal;
}

B1BOptrFDTrackData::~B1BOptrFDTrackData()
{
//	//TODO: fix
//	if ( (fFDState == FDState::toBeSplitCompt) || (fFDState == FDState::toBeSplitRayl) )
//	    {
//	      G4ExceptionDescription ed;
//	      ed << "Track deleted while under G4BOptrForceCollision biasing scheme of operator `";
//	      if ( fFDOperator == nullptr ) ed << "(none)"; else ed << fFDOperator->GetName();
//	      ed <<"'. Will result in inconsistencies.";
//	      G4Exception(" B1BOptrFDTrackData::~B1BOptrFDTrackData()",
//			  "BIAS.GEN.19",
//			  JustWarning,
//			  ed);
//	    }
}

void B1BOptrFDTrackData::Print() const
{
  G4cout << " B1BOptrFDTrackData object : " << this << G4endl;
  G4cout << "     Local Estimation operator : "; if ( fFDOperator == nullptr ) G4cout << "(none)"; else G4cout << fFDOperator->GetName(); G4cout << G4endl;
  G4cout << "     Local Estimation state    : ";
  switch ( fFDState )
    {
    case FDState::normal :
      G4cout << "normal ";
      break;
    case FDState::toBeFreeFlight :
      G4cout << "to be free flight forced ";
      break;
    default:
      break;
    }
  G4cout << G4endl;
}
