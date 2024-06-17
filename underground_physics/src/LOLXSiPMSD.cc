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
/// \file src/LOLXPMTSD.cc
/// \brief Implementation of the LOLXPMTSD class
//
//
#include "LOLXSiPMSD.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LOLXSiPMSD::LOLXSiPMSD(G4String name)
  : G4VSensitiveDetector(name),SiPMHitCollection(0),SiPMPositionsX(0)
  ,SiPMPositionsY(0),SiPMPositionsZ(0)
{
  collectionName.insert("sipmHitCollection");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LOLXSiPMSD::~LOLXSiPMSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LOLXSiPMSD::Initialize(G4HCofThisEvent* hitsCE){
  SiPMHitCollection = new LOLXSiPMHitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  //Store collection with event and keep ID
  static G4int hitCID = -1;
  if(hitCID<0){
    hitCID = GetCollectionID(0);
  }
  hitsCE->AddHitsCollection( hitCID, SiPMHitCollection );
}

LOLXSiPMHit* LOLXSiPMSD::createHit(const G4Step* step) {

  const G4DynamicParticle* particle = step->GetTrack()->GetDynamicParticle();
  G4StepPoint* stepPoint = step->GetPostStepPoint();
  const G4TouchableHistory* touchable = reinterpret_cast<const G4TouchableHistory*>(stepPoint->GetTouchable());

  G4int sipmid =  4*touchable->GetReplicaNumber(1) + touchable->GetReplicaNumber(0);
  //G4int sipmid =  4*touchable->GetCopyNumber(1)+touchable->GetCopyNumber(0);
  //printf("touchable->GetReplicaNumber(1) = %d touchable->GetReplicaNumber(0) = %d\n",touchable->GetReplicaNumber(1),touchable->GetVolume(0)->GetCopyNo());
  G4VPhysicalVolume* physVol = touchable->GetVolume(0);
  //printf("volume = %s\n",physVol->GetName().data());
  // printf("create hit for %d %d %d\n",sipmid,SiPMGlobalID[sipmid],SiPMPackageID[sipmid]);
  //Find the correct hit collection
  LOLXSiPMHit* hit = new LOLXSiPMHit(); //so create new hit
  hit->SetSiPMPhysVol(physVol);
  hit->SetSiPMGlobalID(SiPMGlobalID[sipmid]);
  hit->SetSiPMPackageID(SiPMPackageID[sipmid]);
  hit->SetSiPMChipID(SiPMChipID[sipmid]);
  hit->SetHitTime(stepPoint->GetGlobalTime()/ns);
  hit->SetEKin(particle->GetKineticEnergy()/eV);
  hit->SetPdgId(particle->GetPDGcode());
  hit->SetTrackId(step->GetTrack()->GetTrackID());
  hit->SetWorldPosition(stepPoint->GetPosition());
  hit->SetPosition(touchable->GetHistory()->GetTopTransform().TransformPoint(stepPoint->GetPosition())/mm);
  hit->SetMomentum((m/(g*s))*(touchable->GetHistory()->GetTopTransform().TransformPoint(particle->GetMomentum())));
  if(step->GetTrack()->GetCreatorProcess()){
    if(step->GetTrack()->GetCreatorProcess()->GetProcessName()=="Scintillation") hit->SetHitProcess(1);
    else if(step->GetTrack()->GetCreatorProcess()->GetProcessName()=="Cerenkov") hit->SetHitProcess(0);
    else if(step->GetTrack()->GetCreatorProcess()->GetProcessName()=="LOLXOpBoundary"){
      hit->SetHitProcess(5);
      // G4cout << "Package ID = " <<SiPMPackageID[sipmid] << G4endl;
      // G4cout << "Position " << stepPoint->GetPosition() << G4endl;
    }
   else if(step->GetTrack()->GetCreatorProcess()->GetProcessName()=="LOLXOpWLS"){ 
      hit->SetHitProcess(6);
  } else {hit->SetHitProcess(7);}
  //printf("hit on SiPM %d process = %d\n",sipmid,hit->GetHitProcess());
  hit->SetDrawit(true);
 } 
  return hit;
}

//G4bool LOLXSiPMSD::ProcessHitsForced(const G4Step* step, G4TouchableHistory* ROhist)
G4bool LOLXSiPMSD::ProcessHitsForced(const G4Step* step)
{
  //printf("ProcessHitsForced\n");
    if(step->GetTrack()->GetDefinition()
     != G4OpticalPhoton::OpticalPhotonDefinition()) return false;

    SiPMHitCollection->insert(createHit(step));
    return true;
}

G4bool LOLXSiPMSD::ProcessHits(G4Step* step, G4TouchableHistory* ROhist) {

  //printf("hit sipm\n");

  //try to put something in to make sure photons are not entering from wrong side.

  if(step->GetTrack()->GetDefinition()
     != G4OpticalPhoton::OpticalPhotonDefinition()){
    //printf("not optical photon\n");
    return false;
  }
  
  //Get Edep
  G4double edep = step->GetTotalEnergyDeposit();

  //Print some stuff
  if(step->GetTrack()->GetCreatorProcess()){
  	// if(step->GetTrack()->GetCreatorProcess()->GetProcessName()=="eXT"){
    //   G4cout << "Stopping ext photon ID = " << step->GetTrack()->GetTrackID()<< G4endl;
    // }
    G4StepPoint* pPostStepPoint = step->GetPostStepPoint();
    G4VPhysicalVolume* thePostPV = pPostStepPoint->GetPhysicalVolume();
    G4String postPV_name = thePostPV->GetName();
    G4StepPoint* pPreStepPoint = step->GetPreStepPoint();
    G4VPhysicalVolume* thePrePV = pPreStepPoint->GetPhysicalVolume();
    G4String prePV_name = thePrePV->GetName();
    G4double edep = step->GetTotalEnergyDeposit();

    // G4cout << "creator process = "<< step->GetTrack()->GetCreatorProcess()->GetProcessName()<<G4endl;
    // G4cout << "sd step status = "<< pPostStepPoint->GetStepStatus()<<G4endl;
    // G4cout << "Post PV name = "<< postPV_name <<G4endl;
    // G4cout << "Pre PV name = "<< prePV_name <<G4endl;
    // G4cout << "Step EDEP = "<< edep <<G4endl;

    // // If we don't have edep, we don't have detection!!
    //Return without creating a hit
  }

  if((fabs(edep) - 1e-12) < 0){
    // G4cout << "No EDEP, skipping" <<G4endl;
    // G4cout << "Momentum:"<< step->GetTrack()->GetMomentumDirection()<<G4endl;
    // G4cout << "IsFirstStepInFlag:"<< step->IsFirstStepInVolume() <<G4endl;
    if(!step->IsFirstStepInVolume()) step->GetTrack()->SetTrackStatus(fStopAndKill);
    return true;
  }
  // G4cout << "Creating hit" <<G4endl;
  SiPMHitCollection->insert(createHit(step));

  // Stop track.
  //aParticleChange.ProposeTrackStatus(fSuspend);
  step->GetTrack()->SetTrackStatus(fStopAndKill);
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LOLXSiPMSD::EndOfEvent(G4HCofThisEvent* ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LOLXSiPMSD::clear() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LOLXSiPMSD::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LOLXSiPMSD::PrintAll() {}
