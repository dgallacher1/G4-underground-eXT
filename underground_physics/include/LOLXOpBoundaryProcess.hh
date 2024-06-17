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
// $Id: G4OpBoundaryProcess.hh 104458 2017-05-31 15:54:11Z gcosmo $
//
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Boundary Process Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpBoundaryProcess.hh
// Description: Discrete Process -- reflection/refraction at
//                                  optical interfaces
// Version:     1.1
// Created:     1997-06-18
// Modified:    2005-07-28 add G4ProcessType to constructor
//              1999-10-29 add method and class descriptors
//              1999-10-10 - Fill NewMomentum/NewPolarization in
//                           DoAbsorption. These members need to be
//                           filled since DoIt calls
//                           aParticleChange.SetMomentumChange etc.
//                           upon return (thanks to: Clark McGrew)
//              2006-11-04 - add capability of calculating the reflectivity
//                           off a metal surface by way of a complex index
//                           of refraction - Thanks to Sehwook Lee and John
//                           Hauptman (Dept. of Physics - Iowa State Univ.)
//              2009-11-10 - add capability of simulating surface reflections
//                           with Look-Up-Tables (LUT) containing measured
//                           optical reflectance for a variety of surface
//                           treatments - Thanks to Martin Janecek and
//                           William Moses (Lawrence Berkeley National Lab.)
//              2013-06-01 - add the capability of simulating the transmission
//                           of a dichronic filter
//              2017-02-24 - add capability of simulating surface reflections
//                           with Look-Up-Tables (LUT) developed in DAVIS
//
// Author:      Peter Gumplinger
//              adopted from work by Werner Keil - April 2/96
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef LOLXOpBoundaryProcess_h
#define LOLXOpBoundaryProcess_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "geomdefs.hh"
#include "Randomize.hh"

#include "G4RandomTools.hh"
#include "G4RandomDirection.hh"

#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "LOLXOpticalSurface.hh"
#include "G4OpticalPhoton.hh"
#include "G4TransportationManager.hh"

//LolX members
#include "LOLXExtendedOpticalPhysics.hh"
#include "ExternalCrossTalkProcess.hh"

//ROOT members
#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "Math/Transform3D.h"
#include "TRandom.h"
#include "TMath.h"


// Class Description:
// Discrete Process -- reflection/refraction at optical interfaces.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////
//Put enum into namespace LOLXOpBoundary to prevent collision with G4OpBoundaryProcess
enum LOLXOpBoundaryProcessStatus {  lUndefined,
                                  lTransmission, lFresnelRefraction,
                                  lFresnelReflection, lTotalInternalReflection,
                                  lLambertianReflection, lLobeReflection,
                                  lSpikeReflection, lBackScattering,
                                  lAbsorption, lDetection, lNotAtBoundary,
                                  lSameMaterial, lStepTooSmall, lNoRINDEX,
                                  lPolishedLumirrorAirReflection,
                                  lPolishedLumirrorGlueReflection,
                                  lPolishedAirReflection,
                                  lPolishedTeflonAirReflection,
                                  lPolishedTiOAirReflection,
                                  lPolishedTyvekAirReflection,
                                  lPolishedVM2000AirReflection,
                                  lPolishedVM2000GlueReflection,
                                  lEtchedLumirrorAirReflection,
                                  lEtchedLumirrorGlueReflection,
                                  lEtchedAirReflection,
                                  lEtchedTeflonAirReflection,
                                  lEtchedTiOAirReflection,
                                  lEtchedTyvekAirReflection,
                                  lEtchedVM2000AirReflection,
                                  lEtchedVM2000GlueReflection,
                                  lGroundLumirrorAirReflection,
                                  lGroundLumirrorGlueReflection,
                                  lGroundAirReflection,
                                  lGroundTeflonAirReflection,
                                  lGroundTiOAirReflection,
                                  lGroundTyvekAirReflection,
                                  lGroundVM2000AirReflection,
                                  lGroundVM2000GlueReflection,
                                  lDichroic };
class LOLXOpBoundaryProcess : public G4VDiscreteProcess
{

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        LOLXOpBoundaryProcess(const G4String& processName = "LOLXOpBoundary",
                                     G4ProcessType type = fOptical);
        ~LOLXOpBoundaryProcess();

private:

        LOLXOpBoundaryProcess(const LOLXOpBoundaryProcess &right);

        //////////////
        // Operators
        //////////////

        LOLXOpBoundaryProcess& operator=(const LOLXOpBoundaryProcess &right);

public:

	////////////
        // Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an optical photon.

        G4double GetMeanFreePath(const G4Track& ,
                                 G4double ,
                                 G4ForceCondition* condition);
        // Returns infinity; i. e. the process does not limit the step,
        // but sets the 'Forced' condition for the DoIt to be invoked at
        // every step. However, only at a boundary will any action be
        // taken.

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                        const G4Step&  aStep);
        // This is the method implementing boundary processes.

        LOLXOpBoundaryProcessStatus GetStatus() const;
        // Returns the current status.

        void SetInvokeSD(G4bool );
        // Set flag for call to InvokeSD method.

private:

        G4bool G4BooleanRand(const G4double prob) const;

        G4ThreeVector GetFacetNormal(const G4ThreeVector& Momentum,
                                     const G4ThreeVector&  Normal) const;

        void DielectricMetal();
        void DielectricDielectric();

        void DielectricLUT();
        void DielectricLUTDAVIS();

        void DielectricDichroic();

        void ChooseReflection();
        void DoAbsorption();
        void DoReflection();

        G4double GetIncidentAngle();
        // Returns the incident angle of optical photon

        G4double GetReflectivity(G4double E1_perp,
                                 G4double E1_parl,
                                 G4double incidentangle,
                                 G4double RealRindex,
                                 G4double ImaginaryRindex);
        // Returns the Reflectivity on a metalic surface

        void CalculateReflectivity(void);

        void BoundaryProcessVerbose(void) const;

        // Invoke SD for post step point if the photon is 'detected'
        G4bool InvokeSD(const G4Step* step);
        
private:

        G4double thePhotonMomentum;

        G4ThreeVector OldMomentum;
        G4ThreeVector OldPolarization;

        G4ThreeVector NewMomentum;
        G4ThreeVector NewPolarization;

        G4ThreeVector theGlobalNormal;
        G4ThreeVector theFacetNormal;

        G4Material* Material1;
        G4Material* Material2;

        LOLXOpticalSurface* OpticalSurface;

        G4MaterialPropertyVector* PropertyPointer;
        G4MaterialPropertyVector* PropertyPointer1;
        G4MaterialPropertyVector* PropertyPointer2;

        G4double Rindex1;
        G4double Rindex2;

        G4double cost1, cost2, sint1, sint2;

        LOLXOpBoundaryProcessStatus theStatus;

        LOLXOpticalSurfaceModel theModel;

        LOLXOpticalSurfaceFinish theFinish;

        G4double theReflectivity;
        G4double theEfficiency;
        G4double theTransmittance;
        G4double theSurfaceRoughness;

       //SiPM incident angle efficiency correction
        G4double efficiencyCosThetaCorrection;
        G4double efficiencyCosThetaOffset;


        G4double prob_sl, prob_ss, prob_bs;

        G4int iTE, iTM;

        G4double kCarTolerance;

        size_t idx, idy;
        G4Physics2DVector* DichroicVector;


        size_t idxA, idyA;
        G4Physics2DVector* DichroicVectorA;

        G4bool fInvokeSD;
};

////////////////////
// Inline methods
////////////////////

inline
G4bool LOLXOpBoundaryProcess::G4BooleanRand(const G4double prob) const
{
  /* Returns a random boolean variable with the specified probability */

  return (G4UniformRand() < prob);
}

inline
G4bool LOLXOpBoundaryProcess::IsApplicable(const G4ParticleDefinition&
                                                       aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

inline
LOLXOpBoundaryProcessStatus LOLXOpBoundaryProcess::GetStatus() const
{
   return theStatus;
}

inline
void LOLXOpBoundaryProcess::SetInvokeSD(G4bool flag)
{
  fInvokeSD = flag;
}

inline
void LOLXOpBoundaryProcess::ChooseReflection()
{
                 G4double rand = G4UniformRand();
                 if ( rand >= 0.0 && rand < prob_ss ) {
                    theStatus = lSpikeReflection;
                    theFacetNormal = theGlobalNormal;
                 }
                 else if ( rand >= prob_ss &&
                           rand <= prob_ss+prob_sl) {
                    theStatus = lLobeReflection;
                 }
                 else if ( rand > prob_ss+prob_sl &&
                           rand < prob_ss+prob_sl+prob_bs ) {
                    theStatus = lBackScattering;
                 }
                 else {
                    theStatus = lLambertianReflection;
                 }
}

inline
void LOLXOpBoundaryProcess::DoAbsorption()
{
              theStatus = lAbsorption;
              //Get angle of incidence
              G4double incidentAngle = GetIncidentAngle();
              incidentAngle = -1*(incidentAngle-CLHEP::pi);//Need angle WRT to surface normal
              if(incidentAngle > CLHEP::pi/2) incidentAngle = incidentAngle - CLHEP::pi/2; //Photon is going opposite to surface normal, so flip incidence to other side
              G4double corrEff = efficiencyCosThetaCorrection*cos(incidentAngle)+efficiencyCosThetaOffset;
              corrEff *= theEfficiency;
              
            //   G4cout << Form("Absorbing photon (eff,wl) = (%f,%f)\n",theEfficiency, (CLHEP::h_Planck*CLHEP::c_light/thePhotonMomentum/CLHEP::nm));
              
              if (corrEff <= 0) corrEff = 0.0;
              //if ( G4BooleanRand(corrEff) ) {//Disable for now - Unity Efficiency for all angles and WLs
              if ( G4BooleanRand(theEfficiency) ) {
               //   G4cout << "Detected photon!"<<G4endl;
                 // EnergyDeposited =/= 0 means: photon has been detected
                 theStatus = lDetection;
                 aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
              }
              else {
                 aParticleChange.ProposeLocalEnergyDeposit(0.0);
                 aParticleChange.ProposeTrackStatus(fStopAndKill);//Kill absorbed photons
              }

              NewMomentum = OldMomentum;
              NewPolarization = OldPolarization;
}

inline
void LOLXOpBoundaryProcess::DoReflection()
{        
      // G4cout << "Reflection = " << theStatus <<G4endl;
        if ( theStatus == lLambertianReflection ) {
          NewMomentum = G4LambertianRand(theGlobalNormal);
          theFacetNormal = (NewMomentum - OldMomentum).unit();
        } else if ( theFinish == lground ) {
          theStatus = lLobeReflection;
          if ( PropertyPointer1 && PropertyPointer2 ){
          } else {
             theFacetNormal =
                 GetFacetNormal(OldMomentum,theGlobalNormal);
          }
          G4double PdotN = OldMomentum * theFacetNormal;
          NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
        } else {
         //  G4cout << "Reflecting, old mom = " << Form("[%f,%f,%f]",OldMomentum.x(),OldMomentum.y(),OldMomentum.z()) <<G4endl;
          theStatus = lSpikeReflection;
          theFacetNormal = theGlobalNormal;
          G4double PdotN = OldMomentum * theFacetNormal;
          NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
         //  G4cout << "Reflecting, new mom = " << Form("[%f,%f,%f]",NewMomentum.x(),NewMomentum.y(),NewMomentum.z()) <<G4endl;
        }
        G4double EdotN = OldPolarization * theFacetNormal;
        NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
}


#endif /* LOLXOpBoundaryProcess_h */
