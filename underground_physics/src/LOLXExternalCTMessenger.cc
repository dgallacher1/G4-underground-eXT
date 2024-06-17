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
/// \file src/LOLXExternalCTMessenger.cc
/// \brief Implementation of the LOLXExternalCTMessenger class
//
//
#include "LOLXExternalCTMessenger.hh"

#include "ExternalCrossTalkProcess.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LOLXExternalCTMessenger::LOLXExternalCTMessenger(ExternalCrossTalkProcess* ext)
 : fExternalCrosstalk(ext)
{
  fPhysDir = new G4UIdirectory("/LOLX/eXT/");
  fPhysDir->SetGuidance("Detector Physics Lists Properties Control");

  fVerboseCmd = new G4UIcmdWithAnInteger("/LOLX/eXT/verbose",this);
  fVerboseCmd->SetGuidance("Set the verbosity of the ext processor list.");
  fVerboseCmd->SetParameterName("verbose",true);
  fVerboseCmd->SetDefaultValue(1);

  fSetEXTProbCmd = new G4UIcmdWithADouble("/LOLX/eXT/prob",this);
  fSetEXTProbCmd->SetGuidance("Set the probability for the external cross-talk process");
  fSetEXTProbCmd->SetDefaultValue(1.0);

  fSetEXTIntensityCmd = new G4UIcmdWithADouble("/LOLX/eXT/intensity",this);
  fSetEXTIntensityCmd->SetGuidance("Set the intensity in mean num photons for the external cross-talk process");
  fSetEXTIntensityCmd->SetDefaultValue(0.0);

  fSetEXTSmearCmd = new G4UIcmdWithADouble("/LOLX/eXT/smear_level",this);
  fSetEXTSmearCmd->SetGuidance("Set the amount of gaussian smearing of the mean num photons for the external cross-talk process");
  fSetEXTSmearCmd->SetDefaultValue(0.0);

  fSetEXTTargetNumberCmd = new G4UIcmdWithAnInteger("/LOLX/eXT/target_number",this);
  fSetEXTTargetNumberCmd->SetGuidance("Set the number of the package to trigger the eXT process for, -1 for all");
  fSetEXTTargetNumberCmd->SetDefaultValue(-1);

  fSetEXTDisabledNumberCmd = new G4UIcmdWithAnInteger("/LOLX/eXT/disabled_package",this);
  fSetEXTDisabledNumberCmd->SetGuidance("Set the number of the package to disable the eXT process for, -1 for none");
  fSetEXTDisabledNumberCmd->SetDefaultValue(-1);

  fRefractionBoolCmd= new G4UIcmdWithAnInteger("/LOLX/eXT/noRefraction",this);
  fRefractionBoolCmd->SetGuidance("Set whether to turn off refraction to angular distribution in external cross-talk");
  fRefractionBoolCmd->SetParameterName("noSnells",true);
  fRefractionBoolCmd->SetDefaultValue(0);

  fMethodBoolCmd= new G4UIcmdWithAnInteger("/LOLX/eXT/method",this);
  fMethodBoolCmd->SetGuidance("Set whether to use Kurtis' simulation output or use simple optics");
  fMethodBoolCmd->SetParameterName("method",true);
  fMethodBoolCmd->SetDefaultValue(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LOLXExternalCTMessenger::~LOLXExternalCTMessenger(){
  delete fVerboseCmd;
  delete fSetEXTProbCmd;
  delete fSetEXTIntensityCmd;
  delete fSetEXTTargetNumberCmd;
  delete fRefractionBoolCmd;
  delete fMethodBoolCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LOLXExternalCTMessenger::SetNewValue(G4UIcommand* command, G4String newValue){

  if( command == fVerboseCmd ){
    fExternalCrosstalk->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue));
  }else if( command == fSetEXTProbCmd){
    fExternalCrosstalk->SetEXTProbability(fSetEXTProbCmd->GetNewDoubleValue(newValue));
  }else if( command == fSetEXTIntensityCmd){
    fExternalCrosstalk->SetEXTIntensity(fSetEXTIntensityCmd->GetNewDoubleValue(newValue));
  }else if( command == fSetEXTSmearCmd){
    fExternalCrosstalk->SetEXTSmearing(fSetEXTSmearCmd->GetNewDoubleValue(newValue));
  }else if( command == fSetEXTTargetNumberCmd){
    fExternalCrosstalk->SetTargetPackageNumber(fSetEXTTargetNumberCmd->GetNewIntValue(newValue));
  }else if( command == fSetEXTDisabledNumberCmd){
    fExternalCrosstalk->SetDisabledPackageNumber(fSetEXTDisabledNumberCmd->GetNewIntValue(newValue));
  }else if( command == fRefractionBoolCmd){
    fExternalCrosstalk->SetNoRefractionBool(fRefractionBoolCmd->GetNewIntValue(newValue));
  }else if( command == fMethodBoolCmd){
    G4cout << "LOLXExternalCTMessenger: Set method to "<< newValue <<G4endl;
    fExternalCrosstalk->SetMethodBool(fMethodBoolCmd->GetNewIntValue(newValue));
  }


}
