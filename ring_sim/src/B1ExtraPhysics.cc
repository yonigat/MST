/*
 * B1ExtraPhysics.cc
 *
 *  Created on: Sep 14, 2017
 *      Author: adamgeva
 */
#include "B1ExtraPhysics.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include "G4UserSpecialCuts.hh"
#include "G4StepLimiter.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ExtraPhysics::B1ExtraPhysics()
    : G4VPhysicsConstructor("Extra") { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1ExtraPhysics::~B1ExtraPhysics() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ExtraPhysics::ConstructParticle() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1ExtraPhysics::ConstructProcess()
{
    G4cout << "B1ExtraPhysics:: Add Extra Physics Processes"
              << G4endl;

    auto particleIterator=GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)()) {
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        G4double charge = particle->GetPDGCharge();

        if (!pmanager) {
            std::ostringstream o;
            o << "Particle " << particleName << "without a Process Manager";
            G4Exception("B1ExtraPhysics::ConstructProcess()","",
                         FatalException,o.str().c_str());
        }
        // TODO: break or continue?
        if (particleName == "opticalphoton") break;

        if (charge != 0.0) {
           // All charged particles should have a step limiter
           // to make sure that the steps do not get too long.
           pmanager->AddDiscreteProcess(new G4StepLimiter());
           pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
        }
        if (particleName == "gamma") {
                   // All charged particles should have a step limiter
                   // to make sure that the steps do not get too long.
                   pmanager->AddDiscreteProcess(new G4StepLimiter());
                   pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
        }
//        } else if (particleName == "neutron") {
//          // time cuts for ONLY neutrons:
//          pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
//        } else {
//          // Energy cuts for all other neutral particles
//          pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
//        }
    }
}


