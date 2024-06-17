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
/// \file src/SiPMHit.cc
/// \brief Implementation of the SiPMHit class
//
//
#include "SiPMHit.hh"
#include "G4VTouchable.hh"

G4Allocator<SiPMHit> SiPMHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SiPMHit::SiPMHit(): G4VHit() {hitprocess = -1;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SiPMHit::~SiPMHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SiPMHit::SiPMHit(const SiPMHit &right) : G4VHit()
{
  SiPMGlobalID = right.SiPMGlobalID;
  SiPMPackageID = right.SiPMPackageID;
  SiPMChipID = right.SiPMChipID;
  fPhysVol = right.fPhysVol;
  fDrawit = right.fDrawit;
  hittime = right.hittime;
  hitprocess = right.hitprocess;
  pdgId = right.pdgId;
  trackId = right.trackId;
  eKin = right.eKin;
  position = right.position;
  worldPosition = right.worldPosition;
  momentum = right.momentum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const SiPMHit& SiPMHit::operator=(const SiPMHit &right){

  SiPMGlobalID = right.SiPMGlobalID;
  SiPMPackageID = right.SiPMPackageID;
  SiPMChipID = right.SiPMChipID;
  fPhysVol = right.fPhysVol;
  fDrawit = right.fDrawit;
  hittime = right.hittime;
  hitprocess = right.hitprocess;
  pdgId = right.pdgId;
  trackId = right.trackId;
  eKin = right.eKin;
  position = right.position;
  worldPosition = right.worldPosition;
  momentum = right.momentum; 
  
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SiPMHit::Draw(){
  if(fDrawit&&fPhysVol){ //ReDraw only the PMTs that have hit counts > 0
    //Also need a physical volume to be able to draw anything
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager){//Make sure that the VisManager exists
      G4VisAttributes attribs(G4Colour(1.,0.,0.));
      attribs.SetForceSolid(true);
      G4RotationMatrix rot;
      if(fPhysVol->GetRotation())//If a rotation is defined use it
        rot=*(fPhysVol->GetRotation());
      G4Transform3D trans(rot,fPhysVol->GetTranslation());//Create transform
      pVVisManager->Draw(*fPhysVol,attribs,trans);//Draw it
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SiPMHit::Print() {}
