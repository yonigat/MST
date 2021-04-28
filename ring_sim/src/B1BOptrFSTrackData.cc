//
#include "B1BOptrFSTrackData.hh"
#include "B1BOptrFS.hh"


B1BOptrFSTrackData::B1BOptrFSTrackData( const B1BOptrFS* optr )
: G4VAuxiliaryTrackInformation(),
  fFSOperator( optr ),
  primary(true)
{
  fFSState = FSState::start;
}

B1BOptrFSTrackData::~B1BOptrFSTrackData()
{
//	//TODO: fix
//	if ( (fFSState == FSState::toBeSplitCompt) || (fFSState == FSState::toBeSplitRayl) )
//	    {
//	      G4ExceptionDescription ed;
//	      ed << "Track deleted while under G4BOptrForceCollision biasing scheme of operator `";
//	      if ( fFSOperator == nullptr ) ed << "(none)"; else ed << fFSOperator->GetName();
//	      ed <<"'. Will result in inconsistencies.";
//	      G4Exception(" B1BOptrFSTrackData::~B1BOptrFSTrackData()",
//			  "BIAS.GEN.19",
//			  JustWarning,
//			  ed);
//	    }
}

void B1BOptrFSTrackData::Print() const
{
  G4cout << " B1BOptrFSTrackData object : " << this << G4endl;
  G4cout << "     Local Estimation operator : "; if ( fFSOperator == nullptr ) G4cout << "(none)"; else G4cout << fFSOperator->GetName(); G4cout << G4endl;
  G4cout << "     Local Estimation state    : ";
  switch ( fFSState )
    {
    case FSState::start :
      G4cout << "start biasing ";
      break;
    case FSState::toBeSplitCompt :
      G4cout << "to be cloned ";
      break;
    case FSState::toBeSplitRayl :
      G4cout << "to be cloned ";
      break;
    case FSState::toBeFreeFlight :
      G4cout << "to be free flight forced ";
      break;
    default:
      break;
    }
  G4cout << G4endl;
}
