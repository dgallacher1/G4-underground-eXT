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
/// \file /include/SiPMHit.hh
/// \brief Definition of the SiPMHit class
//
//
#ifndef SiPMHit_h
#define SiPMHit_h 1

#include "MCHit.h"
#include <vector>

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"


class G4VTouchable;

class SiPMHit : public G4VHit
{
  public:
 
    SiPMHit();
    virtual ~SiPMHit();
    SiPMHit(const SiPMHit &right);

    const SiPMHit& operator=(const SiPMHit &right);

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
 
    virtual void Draw();
    virtual void Print();

    inline void SetDrawit(G4bool b){fDrawit=b;}
    inline G4bool GetDrawit(){return fDrawit;}

    inline void SetHitTime(G4double time){hittime = time;}
    G4double GetHitTime(){return hittime;}
    int GetHitProcess(){return hitprocess;}
    inline void SetHitProcess(G4int n){hitprocess = n;}

    inline void SetSiPMGlobalID(G4int n) { SiPMGlobalID = n; }
    inline G4int GetSiPMGlobalID() { return SiPMGlobalID; }

    inline void SetSiPMPackageID(G4int n) { SiPMPackageID = n; }
    inline G4int GetSiPMPackageID() { return SiPMPackageID; }

    inline void SetSiPMChipID(G4int n) { SiPMChipID = n; }
    inline G4int GetSiPMChipID() { return SiPMChipID; }

    inline void SetSiPMPhysVol(G4VPhysicalVolume* physVol){this->fPhysVol=physVol;}
    inline G4VPhysicalVolume* GetSiPMPhysVol(){return fPhysVol;}

    inline void SetPosition(G4ThreeVector _position){position = _position;}
    inline G4ThreeVector GetPosition(){return position;}

    inline void SetWorldPosition(G4ThreeVector _worldPosition){worldPosition = _worldPosition;}
    inline G4ThreeVector GetWorldPosition(){return worldPosition;}

    inline G4ThreeVector GetDirection(){if(momentum.mag()==0.0) return momentum; return (1./momentum.mag())*momentum;}

    inline G4double GetEKin(){return eKin;}
    inline void SetEKin(G4double _eKin){eKin = _eKin;}
    
    inline G4ThreeVector GetMomentum(){return momentum;}
    inline void SetMomentum(G4ThreeVector _momentum){momentum = _momentum;}
    
    inline G4int GetPdgId(){return pdgId;}
    inline void SetPdgId(G4int _pdgId){pdgId = _pdgId;}

    inline G4int GetTrackId(){return trackId;}
    inline void SetTrackId(G4int _trackId){trackId = _trackId;}

    inline void SetMCHit(MCHit* mchit){
        mchit->SetTime(hittime);
        mchit->SetLocalPosition(position.x(),position.y(),position.z());
        mchit->SetLocalMomentum(momentum.x(),momentum.y(),momentum.z());
        mchit->SetProcess(hitprocess);
        mchit->SetSiPMID(SiPMGlobalID);
    }

  private:

    G4int SiPMGlobalID;
    G4int SiPMPackageID;
    G4int SiPMChipID;
    G4VPhysicalVolume* fPhysVol;
    G4bool fDrawit;
    G4double hittime;
    G4int hitprocess;
    G4int pdgId;
    G4int trackId;
    G4double eKin;
    G4ThreeVector position;
    G4ThreeVector worldPosition;
    G4ThreeVector momentum;    
};

typedef G4THitsCollection<SiPMHit> SiPMHitsCollection;

extern G4Allocator<SiPMHit> SiPMHitAllocator;

inline void* SiPMHit::operator new(size_t){
  void *aHit;
  aHit = (void *) SiPMHitAllocator.MallocSingle();
  return aHit;
}

inline void SiPMHit::operator delete(void *aHit){
  SiPMHitAllocator.FreeSingle((SiPMHit*) aHit);
}

#endif
