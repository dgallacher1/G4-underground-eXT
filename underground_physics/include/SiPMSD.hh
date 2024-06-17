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
//
//
#ifndef LOLXSIPMSD_h
#define LOLXSIPMSD_h 1

#include "SiPMHit.hh"
#include "ExtendedOpticalPhysics.hh"

#include <vector>

#include "G4DataVector.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

class G4Step;
class G4HCofThisEvent;

class SiPMSD : public G4VSensitiveDetector
{

  public:

    SiPMSD(G4String name);
    virtual ~SiPMSD();

    virtual void Initialize(G4HCofThisEvent* );
    SiPMHit* createHit(const G4Step* step);
    G4bool ProcessHitsForced(const G4Step* aStep);
    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);

    //A version of processHits that keeps aStep constant
    virtual void EndOfEvent(G4HCofThisEvent* );
    virtual void clear();
    void DrawAll();
    void PrintAll();

    //Initialize the arrays to store sipms positions
    inline void InitSiPMs(G4int nSiPMs){
      if(SiPMPositionsX)delete SiPMPositionsX;
      if(SiPMPositionsY)delete SiPMPositionsY;
      if(SiPMPositionsZ)delete SiPMPositionsZ;
      if(SiPMGlobalID.size()) SiPMGlobalID.clear();
      if(SiPMPackageID.size()) SiPMPackageID.clear();
      if(SiPMChipID.size()) SiPMChipID.clear();
      SiPMPositionsX=new G4DataVector(nSiPMs);
      SiPMPositionsY=new G4DataVector(nSiPMs);
      SiPMPositionsZ=new G4DataVector(nSiPMs);
      SiPMGlobalID= std::vector<int>(nSiPMs,0);
      SiPMPackageID= std::vector<int>(nSiPMs,0);
      SiPMChipID= std::vector<int>(nSiPMs,0);
    }

    //Store a sipm position
    inline void SetSiPMPos(G4int n,G4int package_n,G4int chip_n,G4ThreeVector pos){
      if(SiPMPositionsX)SiPMPositionsX->insertAt(n,pos.x());
      if(SiPMPositionsY)SiPMPositionsY->insertAt(n,pos.y());
      if(SiPMPositionsZ)SiPMPositionsZ->insertAt(n,pos.z());
      if(SiPMGlobalID.size()>(unsigned int)n)SiPMGlobalID[n]=n;
      if(SiPMPackageID.size()>(unsigned int)n)SiPMPackageID[n]=package_n;
      if(SiPMChipID.size()>(unsigned int)n)SiPMChipID[n]=chip_n;
    }

  private:

    SiPMHitsCollection* SiPMHitCollection;

    G4DataVector* SiPMPositionsX;
    G4DataVector* SiPMPositionsY;
    G4DataVector* SiPMPositionsZ;
    std::vector<int> SiPMGlobalID;
    std::vector<int> SiPMPackageID;
    std::vector<int> SiPMChipID;

};

#endif
