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
/// \file include/VUV4Volume.hh
/// \brief Definition of the VUV4Volume class
//
#ifndef VUV4Volume_H
#define VUV4Volume_H 1

#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SystemOfUnits.hh"

#include "TVector3.h"

class VUV4Volume : public G4LogicalVolume
{
  public:

    VUV4Volume();

    G4LogicalVolume* GetChipLogicVolume(){return SiPMlogic;}
    G4ThreeVector GetChipPosition(G4int chipnum,G4RotationMatrix* rotation,G4ThreeVector translation);
    G4double GetWindowFaceZ(){return WindowFaceZ;}
    G4double GetTopFace(){return Package_sizeZ/2.;}

    //For OpticalBoundaries
    std::vector<G4PVPlacement*> SiliconPMs;
    std::vector<G4PVPlacement*> SiO2Layers;
    G4PVPlacement* LXeBuffer;//Buffer layer between SiO2 and QuartzWindow

    void SetSiliconBorder(LOLXOpticalSurface* op);

  private:

    void VisAttributes();
    void SurfaceProperties();
    void GetSiData(std::vector<double>&wl, std::vector<double>&real_index, std::vector<double>&k_index);
    void GetSiO2Data(std::vector<double>&wl, std::vector<double>&real_index, std::vector<double>&k_index);
    void CopyValues();

    G4bool fUpdated;

    G4ThreeVector SiPM_planeG4;
    TVector3 SiPM_plane;
    G4LogicalVolume* SiPMlogic;
    G4LogicalVolume* Ceramiclogic;
    G4LogicalVolume* Windowlogic;
    G4LogicalVolume* SiO2Logic;    
    G4LogicalVolume* BufferLogic;//Buffer volume for optics

    static G4double Package_sizeXY;  // size of the  MPPC package initially made from Xenon.
    static G4double Package_sizeZ;  // everything is housed within this volume
    static G4double Package_border;
    G4double WindowFaceZ;
    G4double SiPM_sizeXY;
    G4double SiPM_sizeZ;

};

#endif
