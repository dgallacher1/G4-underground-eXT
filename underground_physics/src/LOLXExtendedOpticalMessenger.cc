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
/// \file src/LOLXExtendedOpticalMessenger.cc
/// \brief Implementation of the LOLXExtendedOpticalMessenger class
//
//
#include "LOLXExtendedOpticalMessenger.hh"
#include "LOLXOpWLS.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LOLXExtendedOpticalMessenger::LOLXExtendedOpticalMessenger(LOLXOpWLS* event)
 :fOpWLS(event)
{
   fOpDir = new G4UIdirectory("/LOLX/Optical/");
   fOpDir->SetGuidance("Plastic fluoroscence optical parameter control");

   fSetFastWLStime = new G4UIcmdWithADouble("/LOLX/Optical/fastWLStime",this);
   fSetFastWLStime->SetGuidance("Set fast time constant for plastic fluoroscence.");

   fSetSlowWLStime = new G4UIcmdWithADouble("/LOLX/Optical/slowWLStime",this);
   fSetSlowWLStime->SetGuidance("Set slow time constant for plastic fluoroscence.");

   fSetFastWLScomponent = new G4UIcmdWithADouble("/LOLX/Optical/fastWLScomp",this);
   fSetFastWLScomponent->SetGuidance("Set fraction of fast component in fluoroscence.");

   fSetWLSphotonyield = new G4UIcmdWithADouble("LOLX/Optical/WLSphotonyield",this);
   fSetWLSphotonyield->SetGuidance("Set photon yield per absorbed photon in fluoroscence.");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LOLXExtendedOpticalMessenger::~LOLXExtendedOpticalMessenger(){
  delete fOpWLS;
  delete fOpDir;
  delete fSetFastWLStime;
  delete fSetSlowWLStime;
  delete fSetFastWLScomponent;
  delete fSetWLSphotonyield;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LOLXExtendedOpticalMessenger::SetNewValue(G4UIcommand* command, G4String newValue){

    if(command == fSetFastWLStime){
	fOpWLS->SetFastWLStime(fSetFastWLStime->GetNewDoubleValue(newValue));
    }else if(command == fSetSlowWLStime){
	fOpWLS->SetSlowWLStime(fSetSlowWLStime->GetNewDoubleValue(newValue));
    }else if(command == fSetFastWLScomponent){
	fOpWLS->SetFastWLScomp(fSetFastWLScomponent->GetNewDoubleValue(newValue));
    }else if(command == fSetWLSphotonyield){
	fOpWLS->SetWLSphotonyield(fSetWLSphotonyield->GetNewDoubleValue(newValue));
    }

}
