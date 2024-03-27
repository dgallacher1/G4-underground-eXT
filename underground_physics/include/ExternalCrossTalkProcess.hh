/**
 * @file ExternalBoundaryProcess.hh
 * @class ExternalBoundaryProcess
 *
 * @author David Gallacher
 * email: dgallacher@snolab.ca
 * created: 2020-10-07
*/

#ifndef __ExternalCrossTalkProcess__
#define __ExternalCrossTalkProcess__

/////////////
// Includes
/////////////

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4TransportationManager.hh"
#include "G4Alpha.hh"
#include "G4OpticalPhoton.hh"
#include "G4Electron.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"
#include "G4TouchableHistory.hh"
#include "G4NavigationHistory.hh"


#include "LOLXExternalCTMessenger.hh"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "Math/Vector3D.h"
#include "Math/Rotation3D.h"
#include "Math/Transform3D.h"


class LOLXExternalCTMessenger;

class ExternalCrossTalkProcess: public G4VDiscreteProcess
{

 public:
  // Constructor:
  ExternalCrossTalkProcess(const G4String& processName = "eXT",
                      G4ProcessType type = fOptical);
  ~ExternalCrossTalkProcess();


  // Returns true for optical photon and alpha
  G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

  // Returns infinity; i. e. the process does not limit the step,
  // but sets the 'Forced' condition for the DoIt to be invoked at
  // every step. However, only at a boundary will any action be
  // taken.
  G4double GetMeanFreePath(const G4Track&,
                           G4double ,
                           G4ForceCondition* condition);

  // This is the method implementing the process
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack,const G4Step&  aStep);

  //Verbosity
  G4int GetVerboseLevel(void) const;
  void SetVerboseLevel(G4int level);

  //Package ID for target
  G4int GetTargetPackageNumber(){return target_package_number;}
  void SetTargetPackageNumber(G4int _target_package_number){target_package_number=_target_package_number;}

  //Package ID for disabled target
  G4int GetDisabledPackageNumber(){return disabled_package_number;}
  void SetDisabledPackageNumber(G4int _disabled_package_number){disabled_package_number=_disabled_package_number;}


  //Add clipping to emission angle for  the external cross-talk generation?
  int GetNoRefractionBool(){return noSnells;}
  void SetNoRefractionBool(int _noSnells){ noSnells= _noSnells;}
  
  //Decide whether to use simple optics (0) or Kurtis' simulation output
  int GetMethodBool(){return method;}
  void SetMethodBool(int _method){ method = _method;}

  //Static probability for external cross-talk generation
  G4double GetEXTProbability(){return eXTProb;}
  void SetEXTProbability(G4double _eXTProb){eXTProb=_eXTProb;}

  //Mean number of eXT photons to produce
  G4double GetEXTIntensity(){return eXTIntensityMean;}
  void SetEXTIntensity(G4double _eXTIntensityMean){eXTIntensityMean=_eXTIntensityMean;}

  //Smearing of the mean number of eXT photons to produce
  G4double GetEXTSmearing(){return eXTSmearing;}
  void SetEXTSmearing(G4double _eXTSmearing){eXTSmearing = _eXTSmearing;}


  ///Create tracks
  std::vector<G4Track*> CreateEXTPhotons(const G4Step& step);

  ///Random direction aglorithm for 2PI emission
  void GetRandomDir(G4double &x, G4double &y, G4double &z);

  ///Random wavelength 
  G4double GetRandomEnergy();

  ///Add Snells Refraction
  void ApplySnellsLaw(double wl, ROOT::Math::XYZVector &dir, ROOT::Math::XYZVector &norm);

  ///Read ROOT file for WL and refractive index distribution
  void ReadROOTFiles();
  
  ///Apply simulation input from Kurtis'
  void ApplySimInput(double wl, ROOT::Math::XYZVector &dir, ROOT::Math::XYZVector &norm);

  static bool ExternalCrossTalkProcessFlag;


 protected:

  LOLXExternalCTMessenger* messenger;

  G4int target_package_number; /// Package ID number for target, -1 is all packages
  G4int disabled_package_number;/// Package ID for disabled packages, use in combination with -1 for all targets
  int noSnells; ///Option to have no refraction on emission angle for to external cross-talk emission, default false
  G4double eXTProb; /// Static Probability of external cross-talk emission
  G4double eXTIntensityMean;///Mean of a poisson distribution for eXT photon number
  int method; // Option to use Kurtis' simulation output as input if method=1, default=0
  double eXTSmearing;///sigma for gaussian smearing of mean emission, randomly sampled for each event and channel

  TFile *inputFile;//File containing histograms
  TH1D *hEXTWavelength;// TH1D from ROOT file of ext emission histograms 
  TH1D *hSiRefractiveWavelength;// TH1D of Si refractive index vs wavelength 
  TH1D *hSiO2RefractiveWavelength;// TH1D of SiO2 refractive index vs wavelength
  TH1D *hLXeRefractiveWavelength;// TH1D of LXe refractive index vs wavelength

  TFile *inputFileKurtis; // File for input from Kurtis from TRIUMF's EXT simulation data
  TH1D *hEXTWL_K; //TH1D of wavelength to be sampled from 140 nm SiO2 stack file
  TH2D *hTransmission;//Transmission from Si to SiO2
  TH2D *hConversion;//Transmission from SiO2 to LXe, angle refraction

};//End class

////////////////////
// Inline methods
////////////////////
inline
G4bool ExternalCrossTalkProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return ((&aParticleType == G4OpticalPhoton::OpticalPhoton()));
}
inline void ExternalCrossTalkProcess::SetVerboseLevel(int level) { verboseLevel=level; }
inline G4int ExternalCrossTalkProcess::GetVerboseLevel(void) const {return verboseLevel; }


#endif
