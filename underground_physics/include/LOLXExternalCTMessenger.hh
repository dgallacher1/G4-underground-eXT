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
/// Messenger to change properties in LOLXPhysicsList
///
//
//

#ifndef LOLXExternalCTMessenger_h
#define LOLXExternalCTMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

#include "ExternalCrossTalkProcess.hh"

class ExternalCrossTalkProcess;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

class LOLXExternalCTMessenger: public G4UImessenger
{
  public:

    LOLXExternalCTMessenger(ExternalCrossTalkProcess*);
    virtual ~LOLXExternalCTMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:

    ExternalCrossTalkProcess      *fExternalCrosstalk;
    G4UIdirectory        *fPhysDir;
    G4UIcmdWithAnInteger *fVerboseCmd;
    G4UIcmdWithADouble   *fSetEXTProbCmd;
    G4UIcmdWithADouble   *fSetEXTIntensityCmd;
    G4UIcmdWithADouble   *fSetEXTSmearCmd;
    G4UIcmdWithAnInteger *fSetEXTTargetNumberCmd;
    G4UIcmdWithAnInteger *fSetEXTDisabledNumberCmd;
    G4UIcmdWithAnInteger *fRefractionBoolCmd;
    G4UIcmdWithAnInteger *fMethodBoolCmd;
};

#endif
