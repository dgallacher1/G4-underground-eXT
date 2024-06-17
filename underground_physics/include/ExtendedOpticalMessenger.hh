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
/// \file include/ExtendedOpticalMessenger.hh
/// \brief Definition of the ExtendedOpticalMessenger class
/// This class handles information passed from mac file to the optical surface builder
//

#ifndef ExtendedOpticalMessenger_h
#define ExtendedOpticalMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "LOLXOpWLS.hh"

class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIdirectory;
class LOLXOpWLS;

class ExtendedOpticalMessenger : public G4UImessenger
{
  public:

    ExtendedOpticalMessenger(LOLXOpWLS*);
    virtual ~ExtendedOpticalMessenger();

    virtual void SetNewValue(G4UIcommand*,G4String);

  private:

    LOLXOpWLS		*fOpWLS;
    G4UIdirectory	*fOpDir;
    G4UIcmdWithADouble	*fSetFastWLStime;  // To set the fast time constant for plastic fluoroscence
    G4UIcmdWithADouble  *fSetSlowWLStime;  // To set the slow time constant for plastic fluoroscence
    G4UIcmdWithADouble  *fSetFastWLScomponent;  // To set the yield of fast component, slow will be calculated as (1-fast)
    G4UIcmdWithADouble  *fSetWLSphotonyield;  // To set the number of photons emitted per photon absorbed by plastic

};

#endif
