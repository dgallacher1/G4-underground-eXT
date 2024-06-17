//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file src/ExtendedOpticalPhysics.cc
/// \brief Implementation of the ExtendedOpticalPhysics class
//
//
#include "ExtendedOpticalPhysics.hh"

ExtendedOpticalPhysics* ExtendedOpticalPhysics::fExtendedPhysics = NULL;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExtendedOpticalPhysics::ExtendedOpticalPhysics(const G4String& name) : G4VPhysicsConstructor(name)
{
  fExtendedPhysics = this;//Set to this iteration
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExtendedOpticalPhysics::~ExtendedOpticalPhysics() {
  theExternalCrossTalkProcess = NULL;
  verbose = 0;
  fExtendedPhysics = NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Get the singleton to access external crosstalk from boundary process
ExtendedOpticalPhysics* ExtendedOpticalPhysics::GetExtendedOpticalPhysics()
{
  return fExtendedPhysics;
}


void ExtendedOpticalPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExtendedOpticalPhysics::ConstructProcess()
{

  theExternalCrossTalkProcess = new ExternalCrossTalkProcess("eXT",fOptical);
  theExternalCrossTalkProcess->SetVerboseLevel(0);//Verbosity

  fBoundary = new CustomOpBoundaryProcess("CustomOpBoundary",fOptical);//Modified G4OpBoundaryProcess for ExCT

  // Add ExternalCrossTalk Process for SiPMs
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pManager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if(theExternalCrossTalkProcess->IsApplicable(*particle)){
      G4cout << "Adding External Crosstalk process for " << particleName << " \n";
      pManager->AddProcess(theExternalCrossTalkProcess);
      pManager->AddDiscreteProcess(fBoundary);//Modified LOLX G4OpBoundaryProcess
      pManager->SetProcessOrderingToLast(theExternalCrossTalkProcess, idxPostStep);
    }//End of conditional for particle
  }//End of loop over particles
}//end of ConstructProcess()
