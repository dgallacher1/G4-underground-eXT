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
/// \file include/LOLXExtendedOpticalPhysics.hh
/// \brief Definition of the LOLXExtendedOpticalPhysics class
//
// Extended optical physics module to include custom physics processes for opticalphoton's
// Includes - ExternalCrossTalkProcess: SiPM external cross-talk photon generation
//            Produces IR photons from SiPM hits
// Author @ David Gallacher
// Email @ dgallacher@snolab.ca
//
#ifndef LOLXExtendedOpticalPhysics_h
#define LOLXExtendedOpticalPhysics_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4VPhysicsConstructor.hh"
#include <iomanip>
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4OpticalPhoton.hh"


//Custom External cross-talk process
#include "ExternalCrossTalkProcess.hh"
class LOLXOpBoundaryProcess;//Forward declaration to break circular dependency between proceeses
#include "LOLXOpBoundaryProcess.hh"
#include "LOLXOpWLS.hh"


class LOLXExtendedOpticalPhysics : public G4VPhysicsConstructor
{
  private:
    static LOLXExtendedOpticalPhysics* fExtendedPhysics;

  public:
    static LOLXExtendedOpticalPhysics* GetExtendedOpticalPhysics();

  public:

    LOLXExtendedOpticalPhysics(const G4String& name = "extendedOptical");
    virtual ~LOLXExtendedOpticalPhysics();

    // This method will be invoked in the Construct() method.
    // each particle type will be instantiated
    virtual void ConstructParticle();

    // This method will be invoked in the Construct() method.
    // each physics process will be instantiated and
    // registered to the process manager of each particle type
    virtual void ConstructProcess();

    ExternalCrossTalkProcess* GetEXTProcess(){return theExternalCrossTalkProcess;}



  protected:

    ExternalCrossTalkProcess *theExternalCrossTalkProcess;//External Cross-talk process for SiPMs
    LOLXOpBoundaryProcess* fBoundary;//Custom OpBoundaryProcess
    LOLXOpWLS* fLOLXWLS; // LOLX WLS process
    G4int verbose;

};

#endif
