#ifndef B1RegularDetectorConstruction_h
#define B1RegularDetectorConstruction_h 1

#include "globals.hh"
#include "B1DetectorConstruction.hh"

//*******************************************************
/// B1RegularDetectorConstruction class
///
/// Construct the phantom using B1PhantomParameterisatin
///
/// History: 30.11.07  First version
/// \author  P. Arce
//*******************************************************

class B1RegularDetectorConstruction : public B1DetectorConstruction
{
public:

  B1RegularDetectorConstruction();
  ~B1RegularDetectorConstruction();

private:

  virtual void ConstructPhantom();

};

#endif
