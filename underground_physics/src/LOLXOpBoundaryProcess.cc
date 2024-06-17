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
////////////////////////////////////////////////////////////////////////
// Optical Photon Boundary Process Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpBoundaryProcess.cc
// Description: Discrete Process -- reflection/refraction at
//                                  optical interfaces
// Version:     1.1
// Created:     1997-06-18
// Modified:    1998-05-25 - Correct parallel component of polarization
//                           (thanks to: Stefano Magni + Giovanni Pieri)
//              1998-05-28 - NULL Rindex pointer before reuse
//                           (thanks to: Stefano Magni)
//              1998-06-11 - delete *sint1 in oblique reflection
//                           (thanks to: Giovanni Pieri)
//              1998-06-19 - move from GetLocalExitNormal() to the new
//                           method: GetLocalExitNormal(&valid) to get
//                           the surface normal in all cases
//              1998-11-07 - NULL OpticalSurface pointer before use
//                           comparison not sharp for: std::abs(cost1) < 1.0
//                           remove sin1, sin2 in lines 556,567
//                           (thanks to Stefano Magni)
//              1999-10-10 - Accommodate changes done in DoAbsorption by
//                           changing logic in DielectricMetal
//              2001-10-18 - avoid Linux (gcc-2.95.2) warning about variables
//                           might be used uninitialized in this function
//                           moved E2_perp, E2_parl and E2_total out of 'if'
//              2003-11-27 - Modified line 168-9 to reflect changes made to
//                           G4OpticalSurface class ( by Fan Lei)
//              2004-02-02 - Set theStatus = Undefined at start of DoIt
//              2005-07-28 - add G4ProcessType to constructor
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
// 		adopted from work by Werner Keil - April 2/96
// mail:        gum@triumf.ca
//
//
// Modified LoLX Version
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpProcessSubType.hh"

#include "LOLXOpBoundaryProcess.hh"
#include "G4GeometryTolerance.hh"

#include "G4VSensitiveDetector.hh"
#include "G4ParallelWorldProcess.hh"

#include "G4SystemOfUnits.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

// G4OpBoundaryProcess::operator=(const G4OpBoundaryProcess &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

LOLXOpBoundaryProcess::LOLXOpBoundaryProcess(const G4String& processName,
                                               G4ProcessType type)
             : G4VDiscreteProcess(processName, type)
{
        //verboseLevel = 1;
        if ( verboseLevel > 0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }

        SetProcessSubType(fOpBoundary);

        theStatus = lUndefined;
        theModel = lglisur;
        theFinish = lpolished;
        theReflectivity =  1.;
        theEfficiency   =  0.;
        theTransmittance = 0.;
      //  theAbsorb   =  0.;

        theSurfaceRoughness = 0.;

        prob_sl = 0.;
        prob_ss = 0.;
        prob_bs = 0.;

        PropertyPointer  = NULL;
        PropertyPointer1 = NULL;
        PropertyPointer2 = NULL;

        Material1 = NULL;
        Material2 = NULL;

        OpticalSurface = NULL;

        kCarTolerance = G4GeometryTolerance::GetInstance()
                        ->GetSurfaceTolerance();

        iTE = iTM = 0;
        thePhotonMomentum = 0.;
        Rindex1 = Rindex2 = 1.;
        cost1 = cost2 = sint1 = sint2 = 0.;

        idx = idy = 0;
        idxA = idyA = 0;
        DichroicVector = NULL;
        DichroicVectorA = NULL;  // Additional vector for the absorbance
        fInvokeSD = true;

        //Hardcoded SiPM angular dependence fit to data from https://arxiv.org/abs/1910.06438
        efficiencyCosThetaCorrection = 0.529;
        efficiencyCosThetaOffset = 0.4386;
}

// LOLXOpBoundaryProcess::LOLXOpBoundaryProcess(const LOLXOpBoundaryProcess &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

LOLXOpBoundaryProcess::~LOLXOpBoundaryProcess(){}

        ////////////
        // Methods
        ////////////

// PostStepDoIt
// ------------
//

G4VParticleChange*
LOLXOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
        //G4cout << "Start of LOLXOpBoundaryProcess PostStepDoIt"<<G4endl;
        theStatus = lUndefined;
        //verboseLevel=1;
        aParticleChange.Initialize(aTrack);
        aParticleChange.ProposeVelocity(aTrack.GetVelocity());

        // Get hyperStep from  G4ParallelWorldProcess
        //  NOTE: PostSetpDoIt of this process should be
        //        invoked after G4ParallelWorldProcess!

        const G4Step* pStep = &aStep;

        const G4Step* hStep = G4ParallelWorldProcess::GetHyperStep();

        if (hStep) pStep = hStep;

        G4bool isOnBoundary =
                (pStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);

        if (isOnBoundary) {
           Material1 = pStep->GetPreStepPoint()->GetMaterial();
           Material2 = pStep->GetPostStepPoint()->GetMaterial();
        } else {
           theStatus = lNotAtBoundary;
           if ( verboseLevel > 0)
            BoundaryProcessVerbose();
           return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
        }
      //  std::cout<<" Material 1 "<< Material1->GetName() <<" Material 2 "<< Material2->GetName()<<"\n";
        G4VPhysicalVolume* thePrePV  =
                               pStep->GetPreStepPoint() ->GetPhysicalVolume();
        G4VPhysicalVolume* thePostPV =
                               pStep->GetPostStepPoint()->GetPhysicalVolume();

        if ( verboseLevel > 0 )
         {
           G4cout << " Photon at Boundary! " << G4endl;
           if (thePrePV)  G4cout << " thePrePV:  " << thePrePV->GetName()  << G4endl;
           if (thePostPV) G4cout << " thePostPV: " << thePostPV->GetName() << G4endl;
        }

        if (aTrack.GetStepLength()<=kCarTolerance/2){
                theStatus = lStepTooSmall;
                if ( verboseLevel > 0)
                BoundaryProcessVerbose();
                return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
        }

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

        thePhotonMomentum = aParticle->GetTotalMomentum();
        OldMomentum       = aParticle->GetMomentumDirection();
        OldPolarization   = aParticle->GetPolarization();

        if ( verboseLevel > 0 )
        {
           G4cout << " Old Momentum Direction: " << OldMomentum     << G4endl;
           G4cout << " Old Polarization:       " << OldPolarization << G4endl;
        }

        G4ThreeVector theGlobalPoint = pStep->GetPostStepPoint()->GetPosition();

        G4bool valid;
        //  Use the new method for Exit Normal in global coordinates,
        //    which provides the normal more reliably.

        // ID of Navigator which limits step

        G4int hNavId = G4ParallelWorldProcess::GetHypNavigatorID();
        std::vector<G4Navigator*>::iterator iNav =
                G4TransportationManager::GetTransportationManager()->
                                         GetActiveNavigatorsIterator();
        theGlobalNormal =
                   (iNav[hNavId])->GetGlobalExitNormal(theGlobalPoint,&valid);

        if (valid) {
          theGlobalNormal = -theGlobalNormal;
        }
        else
        {
          G4ExceptionDescription ed;
          ed << " LOLXOpBoundaryProcess/PostStepDoIt(): "
                 << " The Navigator reports that it returned an invalid normal"
                 << G4endl;
          G4Exception("LOLXOpBoundaryProcess::PostStepDoIt", "OpBoun01",
                      EventMustBeAborted,ed,
                      "Invalid Surface Normal - Geometry must return valid surface normal");
        }

        if (OldMomentum * theGlobalNormal > 0.0) {
#ifdef G4OPTICAL_DEBUG
           G4ExceptionDescription ed;
           ed << " LOLXOpBoundaryProcess/PostStepDoIt(): "
              << " theGlobalNormal points in a wrong direction. "
              << G4endl;
           ed << "    The momentum of the photon arriving at interface (oldMomentum)"
              << " must exit the volume cross in the step. " << G4endl;
           ed << "  So it MUST have dot < 0 with the normal that Exits the new volume (globalNormal)." << G4endl;
           ed << "  >> The dot product of oldMomentum and global Normal is " << OldMomentum*theGlobalNormal << G4endl;
           ed << "     Old Momentum  (during step)     = " << OldMomentum << G4endl;
           ed << "     Global Normal (Exiting New Vol) = " << theGlobalNormal << G4endl;
           ed << G4endl;
           G4Exception("LOLXOpBoundaryProcess::PostStepDoIt", "OpBoun02",
                       EventMustBeAborted,  // Or JustWarning to see if it happens repeatedbly on one ray
                       ed,
                      "Invalid Surface Normal - Geometry must return valid surface normal pointing in the right direction");
#else
           theGlobalNormal = -theGlobalNormal;
#endif
        }

	G4MaterialPropertiesTable* aMaterialPropertiesTable;
        G4MaterialPropertyVector* Rindex;

	aMaterialPropertiesTable = Material1->GetMaterialPropertiesTable();
        if (aMaterialPropertiesTable) {
		Rindex = aMaterialPropertiesTable->GetProperty(kRINDEX);
	}
	else {
                theStatus = lNoRINDEX;
                if ( verboseLevel > 0)
                 BoundaryProcessVerbose();
                aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                aParticleChange.ProposeTrackStatus(fStopAndKill);
                return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        if (Rindex) {
           Rindex1 = Rindex->Value(thePhotonMomentum);
        }
        else {
	        theStatus = lNoRINDEX;
                if ( verboseLevel > 0)
                 BoundaryProcessVerbose();
                aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                aParticleChange.ProposeTrackStatus(fStopAndKill);
                return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	}

        theReflectivity =  1.;
        theEfficiency   =  0.;
        theTransmittance = 0.;
    //    theAbsorb    = 0.;

        theSurfaceRoughness = 0.;

        theModel = lglisur;
        theFinish = lpolished;

        G4SurfaceType type = dielectric_dielectric;

        Rindex = NULL;
        OpticalSurface = NULL;

        G4LogicalSurface* Surface = NULL;

        Surface = G4LogicalBorderSurface::GetSurface(thePrePV, thePostPV);

        if (Surface == NULL){
          G4bool enteredDaughter= (thePostPV->GetMotherLogical() ==
                                   thePrePV ->GetLogicalVolume());
	  if(enteredDaughter){
	    Surface =
              G4LogicalSkinSurface::GetSurface(thePostPV->GetLogicalVolume());
	    if(Surface == NULL)
	      Surface =
                G4LogicalSkinSurface::GetSurface(thePrePV->GetLogicalVolume());
	  }
	  else {
	    Surface =
              G4LogicalSkinSurface::GetSurface(thePrePV->GetLogicalVolume());
	    if(Surface == NULL)
	      Surface =
                G4LogicalSkinSurface::GetSurface(thePostPV->GetLogicalVolume());
	  }
	}

        if (Surface) OpticalSurface =
           dynamic_cast <LOLXOpticalSurface*> (Surface->GetSurfaceProperty());

        if (OpticalSurface!=NULL) {

           type      = OpticalSurface->GetType();
           theModel  = OpticalSurface->GetModel();
           theFinish = OpticalSurface->GetFinish();

           aMaterialPropertiesTable = OpticalSurface->
                                        GetMaterialPropertiesTable();

           if (aMaterialPropertiesTable) {

              if (theFinish == lpolishedbackpainted ||
                  theFinish == lgroundbackpainted ) {
                  Rindex = aMaterialPropertiesTable->GetProperty(kRINDEX);
	          if (Rindex) {
                     Rindex2 = Rindex->Value(thePhotonMomentum);
                  }
                  else {
                     theStatus = lNoRINDEX;
                     if ( verboseLevel > 0)
                      BoundaryProcessVerbose();
                     aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                     aParticleChange.ProposeTrackStatus(fStopAndKill);
                     return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
                  }
              }

              PropertyPointer =
                      aMaterialPropertiesTable->GetProperty(kREFLECTIVITY);
              PropertyPointer1 =
                      aMaterialPropertiesTable->GetProperty(kREALRINDEX);
              PropertyPointer2 =
                      aMaterialPropertiesTable->GetProperty(kIMAGINARYRINDEX);

              iTE = 1;
              iTM = 1;

              if (PropertyPointer) {

                 theReflectivity =
                          PropertyPointer->Value(thePhotonMomentum);
               //   G4cout << "TheReflectivity is = "<< theReflectivity <<G4endl;

              } else if (PropertyPointer1 && PropertyPointer2) {
                 CalculateReflectivity();

              }

              PropertyPointer =
              aMaterialPropertiesTable->GetProperty(kEFFICIENCY);
              if (PropertyPointer) {
                      theEfficiency =
                      PropertyPointer->Value(thePhotonMomentum);
                     // G4cout << "LOLXOpBoundaryProcess The efficiency = "<< theEfficiency << G4endl;
              }

              PropertyPointer =
              aMaterialPropertiesTable->GetProperty(kTRANSMITTANCE);
              if (PropertyPointer) {
                      theTransmittance =
                      PropertyPointer->Value(thePhotonMomentum);
              }

              if (aMaterialPropertiesTable->
                                     ConstPropertyExists("SURFACEROUGHNESS"))
                 theSurfaceRoughness = aMaterialPropertiesTable->
                                         GetConstProperty(kSURFACEROUGHNESS);

	      if ( theModel == lunified ) {
                 PropertyPointer =
                 aMaterialPropertiesTable->GetProperty(kSPECULARLOBECONSTANT);
                 if (PropertyPointer) {
                         prob_sl =
                         PropertyPointer->Value(thePhotonMomentum);
                 } else {
                         prob_sl = 0.0;
                 }

                 PropertyPointer =
                 aMaterialPropertiesTable->GetProperty(kSPECULARSPIKECONSTANT);
	         if (PropertyPointer) {
                         prob_ss =
                         PropertyPointer->Value(thePhotonMomentum);
                 } else {
                         prob_ss = 0.0;
                 }

                 PropertyPointer =
                 aMaterialPropertiesTable->GetProperty(kBACKSCATTERCONSTANT);
                 if (PropertyPointer) {
                         prob_bs =
                         PropertyPointer->Value(thePhotonMomentum);
                 } else {
                         prob_bs = 0.0;
                 }
              }
           }
           else if (theFinish == lpolishedbackpainted ||
                    theFinish == lgroundbackpainted ) {
                      aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                      aParticleChange.ProposeTrackStatus(fStopAndKill);
                      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
           }
        }

        if (type == dielectric_dielectric ) {
           if (theFinish == lpolished || theFinish == lground ) {
              if (Material1 == Material2){
                 theStatus = lSameMaterial;
                 if ( verboseLevel > 0){
                   BoundaryProcessVerbose();
                 }
		             return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
	            }
              aMaterialPropertiesTable =
                     Material2->GetMaterialPropertiesTable();
              if (aMaterialPropertiesTable)
                 Rindex = aMaterialPropertiesTable->GetProperty(kRINDEX);
              if (Rindex) {
                 Rindex2 = Rindex->Value(thePhotonMomentum);
              }
              else {
                 theStatus = lNoRINDEX;
                 if ( verboseLevel > 0)
                  BoundaryProcessVerbose();
                 aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);
                 aParticleChange.ProposeTrackStatus(fStopAndKill);
                 return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
              }
           }
        }

	if (type == dielectric_metal) {

          DielectricMetal();

	}
        else if (type == dielectric_LUT) {

          DielectricLUT();

        }
        else if (type == dielectric_LUTDAVIS) {

          DielectricLUTDAVIS();

        }
        else if (type == dielectric_dichroic) {

          DielectricDichroic();

        }
        else if (type == dielectric_dielectric) {

          if ( theFinish == lpolishedbackpainted ||
               theFinish == lgroundbackpainted ) {
             DielectricDielectric();
          }
          else {
             G4double rand = G4UniformRand();
            //  G4cout << "theReflectivity ="<<theReflectivity<<G4endl;
             if ( rand > theReflectivity ) {
                if (rand > theReflectivity + theTransmittance) {
                   DoAbsorption();
                } else {
                  // G4cout << "Transmitted 525"<<G4endl;
                   theStatus = lTransmission;
                   NewMomentum = OldMomentum;
                   NewPolarization = OldPolarization;
                }
             }
             else {
               //  G4cout << "Passed Reflectivity Throw"<<G4endl;
                if ( theFinish == lpolishedfrontpainted ) {
                   DoReflection();
                }
                else if ( theFinish == lgroundfrontpainted ) {
                   theStatus = lLambertianReflection;
                   DoReflection();
                }
                else {
                   DielectricDielectric();
                }
             }
          }
        }
        else {

          G4cerr << " Error: G4BoundaryProcess: illegal boundary type " << G4endl;
          return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

        }
        
      //   G4cout << "New Direction! Status = "<< theStatus <<G4endl;
        NewMomentum = NewMomentum.unit();
        NewPolarization = NewPolarization.unit();
      //   G4cout << " New Momentum Direction: " << NewMomentum     << G4endl;

        if ( verboseLevel > 0) {
           G4cout << " New Momentum Direction: " << NewMomentum     << G4endl;
           G4cout << " New Polarization:       " << NewPolarization << G4endl;
           BoundaryProcessVerbose();
        }

        aParticleChange.ProposeMomentumDirection(NewMomentum);
        aParticleChange.ProposePolarization(NewPolarization);

        if ( theStatus == lFresnelRefraction || theStatus == lTransmission ) {
           G4MaterialPropertyVector* groupvel =
           Material2->GetMaterialPropertiesTable()->GetProperty(kGROUPVEL);
           G4double finalVelocity = groupvel->Value(thePhotonMomentum);
           aParticleChange.ProposeVelocity(finalVelocity);
        }

        // G4cout << "Invoke SD flag = "<< fInvokeSD <<G4endl;
        // G4cout << "theStatus flag = ";
        // BoundaryProcessVerbose();
        if ( theStatus == lDetection && fInvokeSD ){
          //Call SDManager to produce hits on SiPMs
          //Call the EXT Process to produce cross-talk photons
         //  G4cout << "Detection!"<<G4endl;
          ExternalCrossTalkProcess* theExtProcess = LOLXExtendedOpticalPhysics::GetExtendedOpticalPhysics()->GetEXTProcess();
          std::vector<G4Track*> tracks = theExtProcess->CreateEXTPhotons(aStep);
          uint nsecondaries = tracks.size();
          for(uint iS = 0; iS < nsecondaries ; iS++){
            aParticleChange.AddSecondary(tracks[iS]);
          }
          InvokeSD(pStep);
          aParticleChange.ProposeTrackStatus(fStopAndKill);//Kill detected photons
        }
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

void LOLXOpBoundaryProcess::BoundaryProcessVerbose() const
{
        if ( theStatus == lUndefined )
                G4cout << " *** Undefined *** " << G4endl;
        if ( theStatus == lTransmission )
                G4cout << " *** Transmission *** " << G4endl;
        if ( theStatus == lFresnelRefraction )
                G4cout << " *** FresnelRefraction *** " << G4endl;
        if ( theStatus == lFresnelReflection )
                G4cout << " *** FresnelReflection *** " << G4endl;
        if ( theStatus == lTotalInternalReflection )
                G4cout << " *** TotalInternalReflection *** " << G4endl;
        if ( theStatus == lLambertianReflection )
                G4cout << " *** LambertianReflection *** " << G4endl;
        if ( theStatus == lLobeReflection )
                G4cout << " *** LobeReflection *** " << G4endl;
        if ( theStatus == lSpikeReflection )
                G4cout << " *** SpikeReflection *** " << G4endl;
        if ( theStatus == lBackScattering )
                G4cout << " *** BackScattering *** " << G4endl;
        if ( theStatus == lPolishedLumirrorAirReflection )
                G4cout << " *** PolishedLumirrorAirReflection *** " << G4endl;
        if ( theStatus == lPolishedLumirrorGlueReflection )
                G4cout << " *** PolishedLumirrorGlueReflection *** " << G4endl;
        if ( theStatus == lPolishedAirReflection )
                G4cout << " *** PolishedAirReflection *** " << G4endl;
        if ( theStatus == lPolishedTeflonAirReflection )
                G4cout << " *** PolishedTeflonAirReflection *** " << G4endl;
        if ( theStatus == lPolishedTiOAirReflection )
                G4cout << " *** PolishedTiOAirReflection *** " << G4endl;
        if ( theStatus == lPolishedTyvekAirReflection )
                G4cout << " *** PolishedTyvekAirReflection *** " << G4endl;
        if ( theStatus == lPolishedVM2000AirReflection )
                G4cout << " *** PolishedVM2000AirReflection *** " << G4endl;
        if ( theStatus == lPolishedVM2000GlueReflection )
                G4cout << " *** PolishedVM2000GlueReflection *** " << G4endl;
        if ( theStatus == lEtchedLumirrorAirReflection )
                G4cout << " *** EtchedLumirrorAirReflection *** " << G4endl;
        if ( theStatus == lEtchedLumirrorGlueReflection )
                G4cout << " *** EtchedLumirrorGlueReflection *** " << G4endl;
        if ( theStatus == lEtchedAirReflection )
                G4cout << " *** EtchedAirReflection *** " << G4endl;
        if ( theStatus == lEtchedTeflonAirReflection )
                G4cout << " *** EtchedTeflonAirReflection *** " << G4endl;
        if ( theStatus == lEtchedTiOAirReflection )
                G4cout << " *** EtchedTiOAirReflection *** " << G4endl;
        if ( theStatus == lEtchedTyvekAirReflection )
                G4cout << " *** EtchedTyvekAirReflection *** " << G4endl;
        if ( theStatus == lEtchedVM2000AirReflection )
                G4cout << " *** EtchedVM2000AirReflection *** " << G4endl;
        if ( theStatus == lEtchedVM2000GlueReflection )
                G4cout << " *** EtchedVM2000GlueReflection *** " << G4endl;
        if ( theStatus == lGroundLumirrorAirReflection )
                G4cout << " *** GroundLumirrorAirReflection *** " << G4endl;
        if ( theStatus == lGroundLumirrorGlueReflection )
                G4cout << " *** GroundLumirrorGlueReflection *** " << G4endl;
        if ( theStatus == lGroundAirReflection )
                G4cout << " *** GroundAirReflection *** " << G4endl;
        if ( theStatus == lGroundTeflonAirReflection )
                G4cout << " *** GroundTeflonAirReflection *** " << G4endl;
        if ( theStatus == lGroundTiOAirReflection )
                G4cout << " *** GroundTiOAirReflection *** " << G4endl;
        if ( theStatus == lGroundTyvekAirReflection )
                G4cout << " *** GroundTyvekAirReflection *** " << G4endl;
        if ( theStatus == lGroundVM2000AirReflection )
                G4cout << " *** GroundVM2000AirReflection *** " << G4endl;
        if ( theStatus == lGroundVM2000GlueReflection )
                G4cout << " *** GroundVM2000GlueReflection *** " << G4endl;
        if ( theStatus == lAbsorption )
                G4cout << " *** Absorption *** " << G4endl;
        if ( theStatus == lDetection )
                G4cout << " *** Detection *** " << G4endl;
        if ( theStatus == lNotAtBoundary )
                G4cout << " *** NotAtBoundary *** " << G4endl;
        if ( theStatus == lSameMaterial )
                G4cout << " *** SameMaterial *** " << G4endl;
        if ( theStatus == lStepTooSmall )
                G4cout << " *** StepTooSmall *** " << G4endl;
        if ( theStatus == lNoRINDEX )
                G4cout << " *** NoRINDEX *** " << G4endl;
        if ( theStatus == lDichroic )
                G4cout << " *** Dichroic Transmission *** " << G4endl;
}

G4ThreeVector
LOLXOpBoundaryProcess::GetFacetNormal(const G4ThreeVector& Momentum,
			            const G4ThreeVector&  Normal ) const
{
        G4ThreeVector FacetNormal;

        if (theModel == lunified || theModel == lLUT || theModel== lDAVIS) {

           /* This function code alpha to a random value taken from the
           distribution p(alpha) = g(alpha; 0, sigma_alpha)*std::sin(alpha),
           for alpha > 0 and alpha < 90, where g(alpha; 0, sigma_alpha)
           is a gaussian distribution with mean 0 and standard deviation
           sigma_alpha.  */

           G4double alpha;

           G4double sigma_alpha = 0.0;
           if (OpticalSurface) sigma_alpha = OpticalSurface->GetSigmaAlpha();

           if (sigma_alpha == 0.0) return FacetNormal = Normal;

           G4double f_max = std::min(1.0,4.*sigma_alpha);

           G4double phi, SinAlpha, CosAlpha, SinPhi, CosPhi, unit_x, unit_y, unit_z;
           G4ThreeVector tmpNormal;

           do {
              do {
                 alpha = G4RandGauss::shoot(0.0,sigma_alpha);
                 // Loop checking, 13-Aug-2015, Peter Gumplinger
              } while (G4UniformRand()*f_max > std::sin(alpha) || alpha >= halfpi );

              phi = G4UniformRand()*twopi;

              SinAlpha = std::sin(alpha);
              CosAlpha = std::cos(alpha);
              SinPhi = std::sin(phi);
              CosPhi = std::cos(phi);

              unit_x = SinAlpha * CosPhi;
              unit_y = SinAlpha * SinPhi;
              unit_z = CosAlpha;

              FacetNormal.setX(unit_x);
              FacetNormal.setY(unit_y);
              FacetNormal.setZ(unit_z);

              tmpNormal = Normal;

              FacetNormal.rotateUz(tmpNormal);
              // Loop checking, 13-Aug-2015, Peter Gumplinger
           } while (Momentum * FacetNormal >= 0.0);
	} else {
           G4double  polish = 1.0;
           if (OpticalSurface) polish = OpticalSurface->GetPolish();

           if (polish < 1.0) {
              do {
                 G4ThreeVector smear;
                 do {
                    smear.setX(2.*G4UniformRand()-1.0);
                    smear.setY(2.*G4UniformRand()-1.0);
                    smear.setZ(2.*G4UniformRand()-1.0);
                    // Loop checking, 13-Aug-2015, Peter Gumplinger
                 } while (smear.mag()>1.0);
                 smear = (1.-polish) * smear;
                 FacetNormal = Normal + smear;
                 // Loop checking, 13-Aug-2015, Peter Gumplinger
              } while (Momentum * FacetNormal >= 0.0);
              FacetNormal = FacetNormal.unit();
           }
           else {
              FacetNormal = Normal;
           }
	}
        return FacetNormal;
}

void LOLXOpBoundaryProcess::DielectricMetal()
{
        G4int n = 0;
        G4double rand, PdotN, EdotN;
        G4ThreeVector A_trans, A_paral;

      //   G4cout << "Here DM: refl = "<< theReflectivity <<G4endl;
      //   G4cout << "TheTransmittance = "<<theTransmittance <<G4endl;

        do {
           n++;
           rand = G4UniformRand();
           if ( rand > theReflectivity && n == 1 ) {
            //   G4cout << "Here 1: R = "<< theReflectivity <<G4endl;
              if (rand > theReflectivity + theTransmittance) {
               //  G4cout << "Here 2: Absorption (Detection)" <<G4endl;
                DoAbsorption();
              } else {
               //  G4cout << "Here 3: Reflection" <<G4endl;
                theStatus = lTransmission;
                NewMomentum = OldMomentum;
                NewPolarization = OldPolarization;
              }
              break;
           } else {
             if (PropertyPointer1 && PropertyPointer2) {
               //  G4cout << "Here 4: Calculate Refl from n,k" <<G4endl;
                if ( n > 1 ) {
                   CalculateReflectivity();
                   if ( !G4BooleanRand(theReflectivity) ) {
                     //  G4cout << "Here 2.1: Absorption (Detection)" <<G4endl;
                      DoAbsorption();
                      break;
                   }
                }
             }

             if ( theModel == lglisur || theFinish == lpolished ) {
               // G4cout << "Here 5: Reflection::Polished" <<G4endl;
               DoReflection();
             } else {

                if ( n == 1 ) ChooseReflection();

                if ( theStatus == lLambertianReflection ) {
                   DoReflection();
                }
                else if ( theStatus == lBackScattering ) {
                   NewMomentum = -OldMomentum;
                   NewPolarization = -OldPolarization;
                }
                else {

                   if(theStatus== lLobeReflection){
                     if ( PropertyPointer1 && PropertyPointer2 ){
                     } else {
                        theFacetNormal =
                            GetFacetNormal(OldMomentum,theGlobalNormal);
                     }
                   }

                   PdotN = OldMomentum * theFacetNormal;
                   NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
                   EdotN = OldPolarization * theFacetNormal;

                   if (sint1 > 0.0 ) {
                      A_trans = OldMomentum.cross(theFacetNormal);
                      A_trans = A_trans.unit();
                   } else {
                      A_trans  = OldPolarization;
                   }
                   A_paral   = NewMomentum.cross(A_trans);
                   A_paral   = A_paral.unit();

                   if(iTE>0&&iTM>0) {
                     NewPolarization =
                           -OldPolarization + (2.*EdotN)*theFacetNormal;
                   } else if (iTE>0) {
                     NewPolarization = -A_trans;
                   } else if (iTM>0) {
                     NewPolarization = -A_paral;
                   }

                }

             }
            //  G4cout << "Here 6: Flip"<<G4endl;
             OldMomentum = NewMomentum;
             OldPolarization = NewPolarization;
	   }
      // G4cout << "Here 7: Around the loop again" <<G4endl; 
      // G4cout << "NewMom*GlobalNormal = "<< NewMomentum * theGlobalNormal <<G4endl;
          // Loop checking, 13-Aug-2015, Peter Gumplinger
	} while (NewMomentum * theGlobalNormal < 0.0);
   // G4cout << "Here 8: End of the road"<<G4endl;
}

void LOLXOpBoundaryProcess::DielectricLUT()
{
        G4int thetaIndex, phiIndex;
        G4double AngularDistributionValue, thetaRad, phiRad, EdotN;
        G4ThreeVector PerpendicularVectorTheta, PerpendicularVectorPhi;

        theStatus = LOLXOpBoundaryProcessStatus(G4int(theFinish) +
                           (G4int(lNoRINDEX)-G4int(lgroundbackpainted)));

        G4int thetaIndexMax = OpticalSurface->GetThetaIndexMax();
        G4int phiIndexMax   = OpticalSurface->GetPhiIndexMax();

        G4double rand;

        do {
           rand = G4UniformRand();
           if ( rand > theReflectivity ) {
              if (rand > theReflectivity + theTransmittance) {
                 DoAbsorption();
              } else {
                 theStatus = lTransmission;
                 NewMomentum = OldMomentum;
                 NewPolarization = OldPolarization;
              }
              break;
           }
           else {
              // Calculate Angle between Normal and Photon Momentum
              G4double anglePhotonToNormal =
                                          OldMomentum.angle(-theGlobalNormal);
              // Round it to closest integer
              G4int angleIncident = G4int(std::floor(180/pi*anglePhotonToNormal+0.5));

              // Take random angles THETA and PHI,
              // and see if below Probability - if not - Redo
              do {
                 thetaIndex = G4RandFlat::shootInt(thetaIndexMax-1);
                 phiIndex = G4RandFlat::shootInt(phiIndexMax-1);
                 // Find probability with the new indeces from LUT
                 AngularDistributionValue = OpticalSurface ->
                   GetAngularDistributionValue(angleIncident,
                                               thetaIndex,
                                               phiIndex);
                // Loop checking, 13-Aug-2015, Peter Gumplinger
              } while ( !G4BooleanRand(AngularDistributionValue) );

              thetaRad = (-90 + 4*thetaIndex)*pi/180;
              phiRad = (-90 + 5*phiIndex)*pi/180;
              // Rotate Photon Momentum in Theta, then in Phi
              NewMomentum = -OldMomentum;

              PerpendicularVectorTheta = NewMomentum.cross(theGlobalNormal);
              if (PerpendicularVectorTheta.mag() < kCarTolerance )
                          PerpendicularVectorTheta = NewMomentum.orthogonal();
              NewMomentum =
                 NewMomentum.rotate(anglePhotonToNormal-thetaRad,
                                    PerpendicularVectorTheta);
              PerpendicularVectorPhi =
                                  PerpendicularVectorTheta.cross(NewMomentum);
              NewMomentum = NewMomentum.rotate(-phiRad,PerpendicularVectorPhi);

              // Rotate Polarization too:
              theFacetNormal = (NewMomentum - OldMomentum).unit();
              EdotN = OldPolarization * theFacetNormal;
              NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
           }
          // Loop checking, 13-Aug-2015, Peter Gumplinger
        } while (NewMomentum * theGlobalNormal <= 0.0);
}

void LOLXOpBoundaryProcess::DielectricLUTDAVIS()
{
  G4int angindex, random, angleIncident;
  G4double ReflectivityValue, elevation, azimuth, EdotN;
  G4double anglePhotonToNormal;

  G4int LUTbin = OpticalSurface->GetLUTbins();

  G4double rand = G4UniformRand();

  do {

     anglePhotonToNormal = OldMomentum.angle(-theGlobalNormal);
     angleIncident = G4int(std::floor(180/pi*anglePhotonToNormal+0.5));

     ReflectivityValue = OpticalSurface -> GetReflectivityLUTValue(angleIncident);

     if ( rand > ReflectivityValue ) {

        if ( theEfficiency > 0 ) {
           DoAbsorption();
           break;
        }
        else {

           theStatus = lTransmission;

           if (angleIncident <= 0.01) {
              NewMomentum = OldMomentum;
              break;

           }

           do {
              random   = G4RandFlat::shootInt(1,LUTbin+1);
              angindex = (((random*2)-1))+angleIncident*LUTbin*2 + 3640000;

              azimuth  = OpticalSurface -> GetAngularDistributionValueLUT(angindex-1);
              elevation= OpticalSurface -> GetAngularDistributionValueLUT(angindex);

           } while ( elevation == 0 && azimuth == 0);

           NewMomentum = -OldMomentum;

           G4ThreeVector v = theGlobalNormal.cross(-NewMomentum);
           G4ThreeVector vNorm = v/v.mag();
           G4ThreeVector u = vNorm.cross(theGlobalNormal);

           u = u *= (std::sin(elevation) * std::cos(azimuth));
           v = vNorm *= (std::sin(elevation) * std::sin(azimuth));
           G4ThreeVector w = theGlobalNormal *= (std::cos(elevation));
           NewMomentum = G4ThreeVector(u+v+w);

           // Rotate Polarization too:
           theFacetNormal = (NewMomentum - OldMomentum).unit();
           EdotN = OldPolarization * theFacetNormal;
           NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
        }
     }
     else {

        theStatus = lLobeReflection;

        if (angleIncident == 0) {
           NewMomentum = -OldMomentum;
           break;
        }

        do {
           random   = G4RandFlat::shootInt(1,LUTbin+1);
           angindex = (((random*2)-1))+(angleIncident-1)*LUTbin*2;

           azimuth   = OpticalSurface -> GetAngularDistributionValueLUT(angindex-1);
           elevation = OpticalSurface -> GetAngularDistributionValueLUT(angindex);
        } while (elevation == 0 && azimuth == 0);

        NewMomentum = -OldMomentum;

        G4ThreeVector v     = theGlobalNormal.cross(-NewMomentum);
        G4ThreeVector vNorm = v/v.mag();
        G4ThreeVector u     = vNorm.cross(theGlobalNormal);

        u = u *= (std::sin(elevation) * std::cos(azimuth));
        v = vNorm *= (std::sin(elevation) * std::sin(azimuth));
        G4ThreeVector w = theGlobalNormal*=(std::cos(elevation));

        NewMomentum = G4ThreeVector(u+v+w);

        // Rotate Polarization too: (needs revision)
        NewPolarization = OldPolarization;
     }
  } while (NewMomentum * theGlobalNormal <= 0.0);
}

void LOLXOpBoundaryProcess::DielectricDichroic()
{
        //G4cout << "In LOLXOpBoundaryProcess::DielectricDichroic" << G4endl;
        // Calculate Angle between Normal and Photon Momentum
        G4double anglePhotonToNormal = OldMomentum.angle(-theGlobalNormal);

        // Round it to closest integer
        G4double angleIncident = std::floor(180/pi*anglePhotonToNormal+0.5);

        //G4cout << "Incident wavelength [nm] = "<<  h_Planck*c_light/thePhotonMomentum/nm << G4endl;

        if (DichroicVector==NULL)
        {
           if (OpticalSurface!=NULL)
           {
            // G4cout << "Reading OpticalSurface (GetDichroicVector)"<< G4endl;
             DichroicVector = OpticalSurface->GetDichroicVector();
             DichroicVectorA = OpticalSurface->GetDichroicVectorA();
           }
        }

        //G4cout << "Incident angle = "<< angleIncident << G4endl;

        if (DichroicVector!=NULL) {
           G4double wavelength = h_Planck*c_light/thePhotonMomentum;
           theTransmittance =
             DichroicVector->Value(wavelength/nm,angleIncident,idx,idy)*perCent;
              // for debug:

		          // G4double absorbance = DichroicVectorA->Value(wavelength/nm,angleIncident,idxA,idyA)*perCent;
              // theAbsorb =DichroicVectorA->Value(wavelength/nm,angleIncident,idxA,idyA)*perCent;
              //G4double theAbsorb= DichroicVectorA->Value(wavelength/nm,angleIncident,idxA,idyA)*perCent;
              // G4cout << " New event***************************************" << G4endl;
              //    G4cout << "wavelength: " << std::floor(wavelength/nm)
              //                             << "nm" << G4endl;
              //    G4cout << "Incident angle: " << angleIncident << "deg" << G4endl;
              //    G4cout << "Transmittance: " << std::floor(theTransmittance/perCent) << "%" << G4endl;
              //    G4cout << " Transmittance is :" << theTransmittance << G4endl;
              //G4double absorb;
          /*  for(int i=100;i<1000;i++){
               wavelength=i;

               G4cout<< " WL="<< wavelength<< G4endl;
               theTransmittance = DichroicVector->Value(wavelength,angleIncident,idx,idy)*perCent;
               G4cout << "Transmittance: " << theTransmittance << G4endl;
               G4double absorb = DichroicVectorA->Value(wavelength,angleIncident,idxA,idyA)*perCent;
                G4cout <<" Absorbance is"<< absorb<<G4endl;

             } */  //<---- THIS was for DEBUG
          //  G4cout << " WL is " << wavelength/nm<< G4endl;
      //      G4cout << "Absorbance is"<< theAbsorb<< G4endl;
        } else {
           G4ExceptionDescription ed;
           ed << " LOLXOpBoundaryProcess/DielectricDichroic(): "
              << " The dichroic surface has no G4Physics2DVector"
              << G4endl;
           G4Exception("LOLXOpBoundaryProcess::DielectricDichroic", "OpBoun03",
                       FatalException,ed,
                       "A dichroic surface must have an associated G4Physics2DVector");
        }

        G4double rand = G4UniformRand();
        // replace the following statement:
        //if ( !G4BooleanRand(theTransmittance) ) { // Not transmitted, so reflect
        // with the one examining rand, since we introduce Absorption

        //G4cout << "Downtown rand = "<< rand << G4endl;

        if( rand> theTransmittance) {
         G4double wavelength= h_Planck*c_light/thePhotonMomentum;
         //G4double absorbance = .6*170*170*170/wavelength/wavelength/wavelength;
         G4double absorbance = 0;
         //G4cout << "Downtown wl = "<< wavelength/nm << G4endl;
         // calculate absorbance  from the DichroicVectorA
         //DichoricVectorA->Dump();
         absorbance = DichroicVectorA->Value(wavelength/nm,angleIncident)*perCent;
          //G4cout << "absorption = " << absorbance << G4endl;
         //G4cout <<" Absorbance is"<< absorbance<< " for WL="<<wavelength/nm<< G4endl;
         //G4cout <<" theTransmittance is"<< theTransmittance<< " for WL="<<wavelength/nm<< G4endl;
          //G4cout <<" theTransmittance + absorbance is"<< theTransmittance+absorbance<< " for WL="<<wavelength/nm<< G4endl;
         //G4cout <<" Rand is "<< rand << G4endl;
         //G4cout << "Downtown absorbance = "<< absorbance << G4endl;
         // absorbance= 0; // uncomment to set absorbance to 0

          if(rand> (theTransmittance+absorbance)) {
            //   G4cout <<" Reflected!! " <<G4endl;
           if ( theModel == lglisur || theFinish == lpolished ) {
              DoReflection();
           } else {
              ChooseReflection();
              if ( theStatus == lLambertianReflection ) {
                 DoReflection();
              } else if ( theStatus == lBackScattering ) {
                 NewMomentum = -OldMomentum;
                 NewPolarization = -OldPolarization;
              } else {
                G4double PdotN, EdotN;
                do {
                   if (theStatus==lLobeReflection)
                      theFacetNormal = GetFacetNormal(OldMomentum,theGlobalNormal);
                   PdotN = OldMomentum * theFacetNormal;
                   NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
                  // Loop checking, 13-Aug-2015, Peter Gumplinger
                } while (NewMomentum * theGlobalNormal <= 0.0);
                EdotN = OldPolarization * theFacetNormal;
                NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;
              }
           }
         }
         else {   aParticleChange.ProposeTrackStatus(fStopAndKill); //DoAbsorption(); //G4cout<< " Absorbed " << G4endl;
       }

       } else {

           theStatus = lDichroic;
           NewMomentum = OldMomentum;
           NewPolarization = OldPolarization;

        }
}

void LOLXOpBoundaryProcess::DielectricDielectric()
{
        G4bool Inside = false;
        G4bool Swap = false;
        G4bool SurfaceRoughnessCriterionPass = 1;
      //   G4cout << "Here DD: transmit = "<< theTransmittance <<G4endl;
      //   G4cout << "Here DD: theSurfaceRoughness = "<< theSurfaceRoughness <<G4endl;
        if (theSurfaceRoughness != 0. && Rindex1 > Rindex2) {
           G4double wavelength = h_Planck*c_light/thePhotonMomentum;
           G4double SurfaceRoughnessCriterion =
             std::exp(-std::pow((4*pi*theSurfaceRoughness*Rindex1*cost1/wavelength),2));
           SurfaceRoughnessCriterionPass =
                                     G4BooleanRand(SurfaceRoughnessCriterion);
        }
           G4double wavelength = h_Planck*c_light/thePhotonMomentum;

        leap:

        G4bool Through = false;
        G4bool Done = false;

        G4double PdotN, EdotN;

        G4ThreeVector A_trans, A_paral, E1pp, E1pl;
        G4double E1_perp, E1_parl;
        G4double s1, s2, E2_perp, E2_parl, E2_total, TransCoeff;
        G4double E2_abs, C_parl, C_perp;
        G4double alpha;

        do {

           if (Through) {
              Swap = !Swap;
              Through = false;
              theGlobalNormal = -theGlobalNormal;
              G4SwapPtr(Material1,Material2);
              G4SwapObj(&Rindex1,&Rindex2);
           }

           if ( theFinish == lpolished ) {
              theFacetNormal = theGlobalNormal;
           }
           else {
              theFacetNormal =
                             GetFacetNormal(OldMomentum,theGlobalNormal);
           }

           PdotN = OldMomentum * theFacetNormal;
           EdotN = OldPolarization * theFacetNormal;

           cost1 = - PdotN;
           if (std::abs(cost1) < 1.0-kCarTolerance){
              sint1 = std::sqrt(1.-cost1*cost1);
              sint2 = sint1*Rindex1/Rindex2;     // *** Snell's Law ***
           }
           else {
              sint1 = 0.0;
              sint2 = 0.0;
           }
           if (sint2 >= 1.0) {

              // Simulate total internal reflection

              if (Swap) Swap = !Swap;

              theStatus = lTotalInternalReflection;

              if ( !SurfaceRoughnessCriterionPass ) theStatus =
                                                       lLambertianReflection;

              if ( theModel == lunified && theFinish != lpolished )
                                                    ChooseReflection();

              if ( theStatus == lLambertianReflection ) {
                 DoReflection();
              }
              else if ( theStatus == lBackScattering ) {
                 NewMomentum = -OldMomentum;
                 NewPolarization = -OldPolarization;
              }
              else {

                 PdotN = OldMomentum * theFacetNormal;
                 NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;
                 EdotN = OldPolarization * theFacetNormal;
                 NewPolarization = -OldPolarization + (2.*EdotN)*theFacetNormal;

              }
           }
           else if (sint2 < 1.0) {

              // Calculate amplitude for transmission (Q = P x N)

              if (cost1 > 0.0) {
                 cost2 =  std::sqrt(1.-sint2*sint2);
              }
              else {
                 cost2 = -std::sqrt(1.-sint2*sint2);
              }

              if (sint1 > 0.0) {
                 A_trans = OldMomentum.cross(theFacetNormal);
                 A_trans = A_trans.unit();
                 E1_perp = OldPolarization * A_trans;
                 E1pp    = E1_perp * A_trans;
                 E1pl    = OldPolarization - E1pp;
                 E1_parl = E1pl.mag();
              }
              else {
                 A_trans  = OldPolarization;
                 // Here we Follow Jackson's conventions and we set the
                 // parallel component = 1 in case of a ray perpendicular
                 // to the surface
                 E1_perp  = 0.0;
                 E1_parl  = 1.0;
              }

              s1 = Rindex1*cost1;
              E2_perp = 2.*s1*E1_perp/(Rindex1*cost1+Rindex2*cost2);
              E2_parl = 2.*s1*E1_parl/(Rindex2*cost1+Rindex1*cost2);
              E2_total = E2_perp*E2_perp + E2_parl*E2_parl;
              s2 = Rindex2*cost2*E2_total;

              if (theTransmittance > 0) TransCoeff = theTransmittance;

              else if (cost1 != 0.0) TransCoeff = s2/s1;
              else TransCoeff = 0.0;

            //   G4cout << "Transmission Coeff = "<<TransCoeff <<G4endl;

              if ( !G4BooleanRand(TransCoeff) ) {
               //   G4cout << "DD Reflection Coeff = "<<1.0-TransCoeff <<G4endl;
               //   G4cout << "Material 1 = " << Material1->GetName()<<G4endl;
               //   G4cout << "Material 2 = " << Material2->GetName()<<G4endl;
               //   // Simulate reflection

                 if (Swap) Swap = !Swap;

                 theStatus = lFresnelReflection;

                 if ( !SurfaceRoughnessCriterionPass ) theStatus =
                                                          lLambertianReflection;

                 if ( theModel == lunified && theFinish != lpolished )
                                                     ChooseReflection();

                 if ( theStatus == lLambertianReflection ) {
                    DoReflection();
                 }
                 else if ( theStatus == lBackScattering ) {
                    NewMomentum = -OldMomentum;
                    NewPolarization = -OldPolarization;
                 }
                 else {

                    PdotN = OldMomentum * theFacetNormal;
                    NewMomentum = OldMomentum - (2.*PdotN)*theFacetNormal;

                    if (sint1 > 0.0) {   // incident ray oblique

                       E2_parl   = Rindex2*E2_parl/Rindex1 - E1_parl;
                       E2_perp   = E2_perp - E1_perp;
                       E2_total  = E2_perp*E2_perp + E2_parl*E2_parl;
                       A_paral   = NewMomentum.cross(A_trans);
                       A_paral   = A_paral.unit();
                       E2_abs    = std::sqrt(E2_total);
                       C_parl    = E2_parl/E2_abs;
                       C_perp    = E2_perp/E2_abs;

                       NewPolarization = C_parl*A_paral + C_perp*A_trans;

                    }

                    else {               // incident ray perpendicular

                       if (Rindex2 > Rindex1) {
                          NewPolarization = - OldPolarization;
                       }
                       else {
                          NewPolarization =   OldPolarization;
                       }

                    }
                 }
              }
              else { // photon gets transmitted
               // G4cout << "Transmitted!"<<G4endl;
                // Simulate transmission/refraction
               // std::cout<<" Transmitted  :  Material 1 "<< Material1->GetName() <<" Material 2 "<< Material2->GetName()<<"\n";
             //   if(Material2->GetName() == "filter_Quartz_BP")
               //     { std::cout<<" Transmittance  :  "<< theTransmittance<<"\n";}
                Inside = !Inside;
                Through = true;
                theStatus = lFresnelRefraction;
                //printf("Transmittance = %f and Photon get transmitted\n",theTransmittance);
                if (sint1 > 0.0) {      // incident ray oblique

                   alpha = cost1 - cost2*(Rindex2/Rindex1);
                   NewMomentum = OldMomentum + alpha*theFacetNormal;
                   NewMomentum = NewMomentum.unit();
//                   PdotN = -cost2;
                   A_paral = NewMomentum.cross(A_trans);
                   A_paral = A_paral.unit();
                   E2_abs     = std::sqrt(E2_total);
                   C_parl     = E2_parl/E2_abs;
                   C_perp     = E2_perp/E2_abs;

                   NewPolarization = C_parl*A_paral + C_perp*A_trans;

                }
                else {                  // incident ray perpendicular

                   NewMomentum = OldMomentum;
                   NewPolarization = OldPolarization;

                }
              }
           }

           OldMomentum = NewMomentum.unit();
           OldPolarization = NewPolarization.unit();

           if (theStatus == lFresnelRefraction) {
              Done = (NewMomentum * theGlobalNormal <= 0.0);
           }
           else {
              Done = (NewMomentum * theGlobalNormal >= -kCarTolerance);
	   }

        // Loop checking, 13-Aug-2015, Peter Gumplinger
	} while (!Done);

        if (Inside && !Swap) {
          if( theFinish == lpolishedbackpainted ||
              theFinish == lgroundbackpainted ) {

              G4double rand = G4UniformRand();
            //   printf("Random = %f\n",rand);
              if ( rand > theReflectivity ) {
                 if (rand > theReflectivity + theTransmittance) {
                    DoAbsorption();
                 } else {
                    theStatus = lTransmission;
                    NewMomentum = OldMomentum;
                    NewPolarization = OldPolarization;
                 }
              }
	      else {
                 if (theStatus != lFresnelRefraction ) {
                    theGlobalNormal = -theGlobalNormal;
                 }
                 else {
                    Swap = !Swap;
                    G4SwapPtr(Material1,Material2);
                    G4SwapObj(&Rindex1,&Rindex2);
                 }
                 if ( theFinish == lgroundbackpainted )
                                        theStatus = lLambertianReflection;

                 DoReflection();

                 theGlobalNormal = -theGlobalNormal;
                 OldMomentum = NewMomentum;

                 goto leap;
              }
          }
        }
}

// GetMeanFreePath
// ---------------
//
G4double LOLXOpBoundaryProcess::GetMeanFreePath(const G4Track& ,
                                              G4double ,
                                              G4ForceCondition* condition)
{
  *condition = Forced;

  return DBL_MAX;
}

G4double LOLXOpBoundaryProcess::GetIncidentAngle()
{
  G4double PdotN = OldMomentum * theFacetNormal;
  G4double magP= OldMomentum.mag();
  G4double magN= theFacetNormal.mag();
  G4double incidentangle = pi - std::acos(PdotN/(magP*magN));

  return incidentangle;
}

G4double LOLXOpBoundaryProcess::GetReflectivity(G4double E1_perp,
                                              G4double E1_parl,
                                              G4double incidentangle,
                                              G4double RealRindex,
                                              G4double ImaginaryRindex)
{
  G4complex Reflectivity, Reflectivity_TE, Reflectivity_TM;
  G4complex N1(Rindex1, 0), N2(RealRindex, ImaginaryRindex);
  G4complex CosPhi;

  G4complex u(1,0);           //unit number 1

  G4complex numeratorTE;      // E1_perp=1 E1_parl=0 -> TE polarization
  G4complex numeratorTM;      // E1_parl=1 E1_perp=0 -> TM polarization
  G4complex denominatorTE, denominatorTM;
  G4complex rTM, rTE;

  G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                    Material1->GetMaterialPropertiesTable();
  G4MaterialPropertyVector* aPropertyPointerR =
                      aMaterialPropertiesTable->GetProperty(kREALRINDEX);
  G4MaterialPropertyVector* aPropertyPointerI =
                      aMaterialPropertiesTable->GetProperty(kIMAGINARYRINDEX);
  if (aPropertyPointerR && aPropertyPointerI) {
     G4double RRindex = aPropertyPointerR->Value(thePhotonMomentum);
     G4double IRindex = aPropertyPointerI->Value(thePhotonMomentum);
     N1 = G4complex(RRindex,IRindex);
  }

  // Following two equations, rTM and rTE, are from: "Introduction To Modern
  // Optics" written by Fowles

  CosPhi=std::sqrt(u-((std::sin(incidentangle)*std::sin(incidentangle))*(N1*N1)/(N2*N2)));

  numeratorTE   = N1*std::cos(incidentangle) - N2*CosPhi;
  denominatorTE = N1*std::cos(incidentangle) + N2*CosPhi;
  rTE = numeratorTE/denominatorTE;

  numeratorTM   = N2*std::cos(incidentangle) - N1*CosPhi;
  denominatorTM = N2*std::cos(incidentangle) + N1*CosPhi;
  rTM = numeratorTM/denominatorTM;

  // This is my calculaton for reflectivity on a metalic surface
  // depending on the fraction of TE and TM polarization
  // when TE polarization, E1_parl=0 and E1_perp=1, R=abs(rTE)^2 and
  // when TM polarization, E1_parl=1 and E1_perp=0, R=abs(rTM)^2

  Reflectivity_TE =  (rTE*conj(rTE))*(E1_perp*E1_perp)
                    / (E1_perp*E1_perp + E1_parl*E1_parl);
  Reflectivity_TM =  (rTM*conj(rTM))*(E1_parl*E1_parl)
                    / (E1_perp*E1_perp + E1_parl*E1_parl);
  Reflectivity    = Reflectivity_TE + Reflectivity_TM;

  do {
     if(G4UniformRand()*real(Reflectivity) > real(Reflectivity_TE))
       {iTE = -1;}else{iTE = 1;}
     if(G4UniformRand()*real(Reflectivity) > real(Reflectivity_TM))
       {iTM = -1;}else{iTM = 1;}
    // Loop checking, 13-Aug-2015, Peter Gumplinger
  } while(iTE<0&&iTM<0);

  return real(Reflectivity);

}

void LOLXOpBoundaryProcess::CalculateReflectivity()
{
  G4double RealRindex =
           PropertyPointer1->Value(thePhotonMomentum);
  G4double ImaginaryRindex =
           PropertyPointer2->Value(thePhotonMomentum);

  // calculate FacetNormal
  if ( theFinish == lground ) {
     theFacetNormal =
               GetFacetNormal(OldMomentum, theGlobalNormal);
  } else {
     theFacetNormal = theGlobalNormal;
  }

  G4double PdotN = OldMomentum * theFacetNormal;
  cost1 = -PdotN;

  if (std::abs(cost1) < 1.0 - kCarTolerance) {
     sint1 = std::sqrt(1. - cost1*cost1);
  } else {
     sint1 = 0.0;
  }

  G4ThreeVector A_trans, A_paral, E1pp, E1pl;
  G4double E1_perp, E1_parl;

  if (sint1 > 0.0 ) {
     A_trans = OldMomentum.cross(theFacetNormal);
     A_trans = A_trans.unit();
     E1_perp = OldPolarization * A_trans;
     E1pp    = E1_perp * A_trans;
     E1pl    = OldPolarization - E1pp;
     E1_parl = E1pl.mag();
  }
  else {
     A_trans  = OldPolarization;
     // Here we Follow Jackson's conventions and we set the
     // parallel component = 1 in case of a ray perpendicular
     // to the surface
     E1_perp  = 0.0;
     E1_parl  = 1.0;
  }

  //calculate incident angle
  G4double incidentangle = GetIncidentAngle();

  //calculate the reflectivity depending on incident angle,
  //polarization and complex refractive

  theReflectivity =
             GetReflectivity(E1_perp, E1_parl, incidentangle,
                                                 RealRindex, ImaginaryRindex);
   // G4cout << Form("CalculateReflectivity: R at %f eV = %f",thePhotonMomentum/1e-6,theReflectivity)<<G4endl;
}

G4bool LOLXOpBoundaryProcess::InvokeSD(const G4Step* pStep)
{
  G4Step aStep = *pStep;

  aStep.AddTotalEnergyDeposit(thePhotonMomentum);

  G4VSensitiveDetector* sd = aStep.GetPostStepPoint()->GetSensitiveDetector();
  if (sd) return sd->Hit(&aStep);
  else return false;
}
