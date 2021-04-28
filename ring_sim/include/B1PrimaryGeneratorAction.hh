#ifndef B1PrimaryGeneratorAction_h
#define B1PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "params.hh"

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <random>
#include <bits/random.h>
#include "ParamsSingleton.hh"


class G4ParticleGun;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
/// The default kinematic is a 80 keV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

class B1PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    B1PrimaryGeneratorAction();    
    virtual ~B1PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);         
  
    // method to access particle gun
    void readSpectrum();
    void readCode();
    G4int getEnergyInd();
  
  private:
    ParamsSingleton* pp = ParamsSingleton::Instance();
    G4ParticleGun* fParticleGun;
//    G4Box* fEnvelopeBox;
    G4int frunIDRand;
    G4int fAllSources;
    G4int** code;
    G4int** intens_code;
    G4int sum_dose;
    G4int num_phot;


//	G4int fspect_test[NUM_OF_SPECTRUM_BINS]; // testing the output spectrum
//	std::ofstream foutput_spect; // file for writing sample spectrum
  protected:
	// create the distribution
	std::vector<G4double> fweights;
	
	std::default_random_engine fgenerator;




};


#endif
