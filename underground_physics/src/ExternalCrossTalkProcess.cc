#include "ExternalCrossTalkProcess.hh"
#include <G4TransportationManager.hh>
#include <G4Navigator.hh>
#include "G4Track.hh"
#include "TString.h"
#include "TMath.h"

#include "TRandom.h"

G4bool ExternalCrossTalkProcess::ExternalCrossTalkProcessFlag = 1;

ExternalCrossTalkProcess::ExternalCrossTalkProcess(
  const G4String& processName,
  G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  //Process init
  //Get the cross-talk probability and other variables

  messenger = new LOLXExternalCTMessenger(this);

  //Defaults, change in macro for other values
  target_package_number = -1; //if -1 use any package
  disabled_package_number = -1; //-1 for no disabled package
  eXTProb = 1.0; //Read in macro
  eXTSmearing = 0.0;//Read in macro
  eXTIntensityMean = 0.0; //Mean of poisson distrbution for photon generation
  noSnells = 0;//default off-> refract, only applies to method=0
  method = 1;//Method = 1 = Kurtis, Method = 0 = Simple optics
  ReadROOTFiles();

}

ExternalCrossTalkProcess::~ExternalCrossTalkProcess()
{
  
  delete gLXeRefractiveWavelength;
  delete gSiRefractiveWavelength;
  delete gSiO2RefractiveWavelength;
  delete hEXTWavelength;
  inputFile->Close();
  delete inputFile;

}

G4VParticleChange*
ExternalCrossTalkProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{

  pParticleChange->Initialize(aTrack);
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();

  //Check if we just passed the boundary
  if (pPostStepPoint->GetStepStatus()) {
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }



  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// GetMeanFreePath
G4double ExternalCrossTalkProcess::GetMeanFreePath(const G4Track& ,
    G4double ,
    G4ForceCondition* condition)
{
  *condition = Forced;

  return DBL_MAX;
}


//Workaround to have photons produced when hits occur
std::vector<G4Track*> ExternalCrossTalkProcess::CreateEXTPhotons(const G4Step& step)
{
  std::vector<G4Track*> ext_photons;
  if(verboseLevel > 0 ) G4cout << "Create EXTPhotons: Prob = "<< eXTProb  << " int = "<< eXTIntensityMean<< " target = "<< target_package_number << " Disabled Package = "<< disabled_package_number<< " noSnells = " << noSnells << " method = "<<method <<G4endl;
  G4Track *track = step.GetTrack();
  G4StepPoint* pPostStepPoint = step.GetPostStepPoint();
  G4StepPoint* pPreStepPoint = step.GetPreStepPoint();
  G4VPhysicalVolume* thePostPV = pPostStepPoint->GetPhysicalVolume();
  G4String postPV_name = thePostPV->GetName();
  G4VPhysicalVolume* thePrePV = pPreStepPoint->GetPhysicalVolume();
  G4String prePV_name = thePrePV->GetName();
  //If the target_package_number is -1, aka default, then we want all packages (default behaviour for when this class is called from LoLXOpBoundary)
  //Otherwise look in the navigation history for the package copy number to match the target package
  if(target_package_number != -1 || disabled_package_number != -1){
    const G4TouchableHistory* touchable_hist = reinterpret_cast<const G4TouchableHistory*>(pPostStepPoint->GetTouchable());
    const G4NavigationHistory* history = touchable_hist->GetHistory();
    G4NavigationHistory *hEditable = new G4NavigationHistory(*history);
    G4int n_steps = 3;//Number of history steps to take, should be in recent history, since this is only called from within the Package in LoLXOpBoundary
    //Loop through nav history to find package number to match target
    for(G4int ih = 0 ; ih < n_steps ; ++ih){
      G4int curr_depth = hEditable->GetDepth();
      G4String volume_name = hEditable->GetVolume(curr_depth)->GetName();
      if(volume_name.contains("Package")){//Check if the current volume in history is a package
        G4int copyNumber = hEditable->GetVolume(curr_depth)->GetCopyNo();//CopyNo
        if(copyNumber == target_package_number){
          if(verboseLevel > 0) G4cout << Form("Matched Target %d to Copy %d, making eXT photons \n",target_package_number,copyNumber);
        }else if(disabled_package_number != -1){
          //Add option to skip if disabled
            if (copyNumber == disabled_package_number){
              if(verboseLevel > 0) G4cout << Form("Disabled Target %d to Copy %d, no photons :( \n",disabled_package_number,copyNumber);
              return ext_photons;//If we've disabled this package, return empty
            } //else do nothing and make photons
        }else{
          return ext_photons;//Return empty, we didn't match IDs
        }
      }
      //Left for debugging
      //G4cout << "Depth = "<< curr_depth <<G4endl;
      //G4cout << Form("Step[%d]: Volume name = ",ih) << hEditable->GetVolume(curr_depth)->GetName() << ", copy# = "<< hEditable->GetVolume(curr_depth)->GetCopyNo()<<G4endl;
      hEditable->BackLevel();//Move back
    }
    delete hEditable;//Cleanup
  }
 

  //If a photon passes into a SiPM then sample for a random number and compare to probability and make new photons
  if ( eXTProb > G4UniformRand() ){
      G4double mean_intensity = eXTIntensityMean;
      //Add smearing if specified
      if(eXTSmearing > 0){
        mean_intensity = gRandom->Gaus(eXTIntensityMean,eXTSmearing*eXTIntensityMean);
      }
      if(verboseLevel > 1) G4cout << Form("eXT Intensity before: %f, and after %f smearing at %f level",eXTIntensityMean,mean_intensity,eXTSmearing)<<G4endl;

      G4double num_ext_photons = gRandom->Poisson(mean_intensity);//Poisson distributed
      //Testing double emission for backwards correlation systematic
      // if(gRandom->Uniform() < 0.1) num_ext_photons *= 2.0;//10% diOCT prob to have double gain
      if (verboseLevel > 0){
        // G4cout << "EXT Post PV name = "<< postPV_name <<G4endl;
        G4cout << "N = "<< num_ext_photons <<" external cross-talk photon generated" <<G4endl;
      }

      //Get the surface normal to the solid we're hitting
      G4VSolid *postSolid = thePostPV->GetLogicalVolume()->GetSolid();
      const G4ThreeVector postPos = track->GetPosition();//Global position
      // G4cout << "Global PostPos = "<< Form("[%f,%f,%f] \n",postPos.x(),postPos.y(),postPos.z());
      //Local postpos
      G4ThreeVector localPos = pPostStepPoint->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(postPos);   
      // G4cout << "Local PostPos = "<< Form("[%f,%f,%f] \n",localPos.x(),localPos.y(),localPos.z());

      G4ThreeVector surfaceNormal = postSolid->SurfaceNormal(postPos);
      // G4cout << "SurfaceNormal = "<< Form("[%f,%f,%f] \n",surfaceNormal.x(),surfaceNormal.y(),surfaceNormal.z());
      G4ThreeVector surfaceNormalLoc = postSolid->SurfaceNormal(localPos);
      // G4cout << "LocalSurfaceNormal = "<< Form("[%f,%f,%f] \n",surfaceNormalLoc.x(),surfaceNormalLoc.y(),surfaceNormalLoc.z());

      //Global direction for local surface normal
      G4ThreeVector globalNormal = pPostStepPoint->GetTouchable()->GetHistory()->GetTopTransform().Inverse().TransformAxis(surfaceNormalLoc).unit();
      // G4cout << "Global SurfaceNormal = "<< Form("[%f,%f,%f] \n",globalNormal.x(),globalNormal.y(),globalNormal.z());

      //Use ROOTs updated vector transformation libraries
      ROOT::Math::XYZVector surfNormal(globalNormal.x(),globalNormal.y(),globalNormal.z());
      surfNormal = surfNormal.Unit();//Normalize vector
      
      //Create a random vector to cross with Surface normal to find a perpendicular vector for X'
      ROOT::Math::XYZVector randVector(G4UniformRand(),G4UniformRand(),G4UniformRand());
      ROOT::Math::XYZVector x_prime = surfNormal.Cross(randVector);
      x_prime = x_prime.Unit();

      //Second perpendicular axis (Y')
      ROOT::Math::XYZVector y_prime = surfNormal.Cross(x_prime);
      y_prime = y_prime.Unit();

      //Rotation matrix to transform our point into the new coordinate system defined by the surface normal
      ROOT::Math::Rotation3D rot3dprime(
      x_prime.X(),x_prime.Y(),x_prime.Z(),
      surfNormal.X(),surfNormal.Y(),surfNormal.Z(),
      y_prime.X(),y_prime.Y(),y_prime.Z());
      rot3dprime = rot3dprime.Inverse();

      //Loop over the number of external cross-talk photons and generate random photons
      // with random directions in a hemisphere about the surface normal
      for(G4int iPh = 0 ; iPh < num_ext_photons ; iPh++){

        G4double photonWL = GetRandomEnergy();// Sample random wl from distribution from TRIUMF
        //Get a random photon in 2pi
        G4double px,py,pz;
        GetRandomDir(px,py,pz);
        ROOT::Math::XYZVector rV(px,py,pz);       

        //Transform to X'Y'Z'
        rV = rot3dprime*rV;
        //Add refraction from snells laws to simulate Si->SiO2->LXe interfaces since we're starting photons in LXe
        if(method == 0){
          ApplySnellsLaw(photonWL,rV,surfNormal);
        }else if(method==1){
          //Use Kurtis's simulation output instead
          // ApplySimInput(photonWL,rV,surfNormal);//Old Bugged method
          ApplySimInputSnells(photonWL,rV,surfNormal);//Updated method using Snells law for refraction, but with transmittance from Kurtis
        }

        //If we have selected it, don't apply refraction
        if(noSnells){
          //Simulate uniform hemisphere for outgoing photons
            GetRandomDir(px,py,pz);
            rV.SetX(px);
            rV.SetY(py);
            rV.SetZ(pz);
            rV = rot3dprime*rV;
        }

        if(rV.Z() < -99.99) {
          // if(verboseLevel>0) G4cout << "Internal reflection, skipping photon" <<G4endl;
          continue;//Catch flagged tracks, set to -100.0 for "reflected" photons
        }else{
          if(verboseLevel>0) G4cout <<"Making EXT Photon"<<G4endl;
        }

        // Create photon momentum direction vector
        G4ParticleMomentum photonMomentum(rV.X(), rV.Y(), rV.Z());

        //Full 4pi for polarization
        G4double phi = 2.0*TMath::Pi()*G4UniformRand();
        G4double sinp = std::sin(phi);
        G4double cosp = std::cos(phi);

        // Determine polarization of new photon
        G4double cost = 1. - 2.0*gRandom->Uniform();
        G4double sint = std::sqrt((1.-cost)*(1.+cost));
        G4double sx = cost*cosp;
        G4double sy = cost*sinp;
        G4double sz = -sint;

        G4ThreeVector photonPolarization(sx, sy, sz);
        G4ThreeVector perp = photonMomentum.cross(photonPolarization);
        phi = 2*TMath::Pi()*G4UniformRand();
        sinp = std::sin(phi);
        cosp = std::cos(phi);
        photonPolarization = cosp * photonPolarization + sinp * perp;
        photonPolarization = photonPolarization.unit();

        //Create cross-talk photon
        G4DynamicParticle *aXTPhoton = new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(),photonMomentum);
        aXTPhoton->SetPolarization(photonPolarization.x(),photonPolarization.y(),photonPolarization.z());
        double photonE_g4 = 2*TMath::Pi()*CLHEP::hbarc/(photonWL*CLHEP::nanometer);//Convert WL to energy in G4 units
        aXTPhoton->SetKineticEnergy(photonE_g4);

        //Shift the produced photons to past the SiO2 window
        G4double window_thickness = 0.20*CLHEP::mm;//110*CLHEP::nm;
        G4ThreeVector positionShift(globalNormal.x(),globalNormal.y(),globalNormal.z());
        positionShift.setMag(1.0*window_thickness);

        //Emission time
        G4double timeE = track->GetGlobalTime();
        G4ThreeVector newPos = track->GetPosition();
        // G4cout << Form("Original Position:[%f,%f,%f] \n",newPos.x(),newPos.y(),newPos.z());
        newPos = newPos + positionShift; //Shift forward 
        // G4cout << Form("Shifted Position:[%f,%f,%f] \n",newPos.x(),newPos.y(),newPos.z()); 
        //New track is daughter of original photon
        G4Track *aSecondaryTrack = new G4Track(aXTPhoton,timeE,newPos);
        aSecondaryTrack->SetParentID(track->GetTrackID());
        aSecondaryTrack->SetTouchableHandle(track->GetTouchableHandle());
        aSecondaryTrack->SetCreatorProcess(this);
        ext_photons.push_back(aSecondaryTrack);
      }//End loop over generated photons
    }//End conditional for external cross-talk photon generation

  return ext_photons;
}


//Simualate random photon direction in 2PI for EXT
void ExternalCrossTalkProcess::GetRandomDir(G4double &x,G4double &y, G4double &z)
{
  //Full 2pi about XYZ
  G4double cost = 1. - 2.0*gRandom->Uniform();
  G4double sint = std::sqrt((1.-cost)*(1.+cost));

  G4double phi = TMath::Pi()*gRandom->Uniform();
  G4double sinp = std::sin(phi);
  G4double cosp = std::cos(phi);

  x = sint*cosp;
  y = sint*sinp;
  z = cost;

}

// Sample random WL from distribuion
G4double ExternalCrossTalkProcess::GetRandomEnergy()
{
  G4double wl;
  if(method==0){
    //Sample random wavelength from distribution. Can use other methods if this is too slow
    wl = hEXTWavelength->GetRandom(); 
    //Edge handling to avoid extrapolations
    if(wl < hEXTWavelength->GetBinLowEdge(1)) wl = hEXTWavelength->GetBinLowEdge(1);
    if(wl > hEXTWavelength->GetBinCenter(hEXTWavelength->GetNbinsX())) wl = hEXTWavelength->GetBinLowEdge(hEXTWavelength->GetNbinsX());
  }else if(method==1){
    wl = hEXTWL_K->GetRandom(); 
    //Edge handling to avoid extrapolations
    if(wl < hEXTWL_K->GetBinLowEdge(1)) wl = hEXTWL_K->GetBinLowEdge(1);
    if(wl > hEXTWL_K->GetBinCenter(hEXTWL_K->GetNbinsX())) wl = hEXTWL_K->GetBinLowEdge(hEXTWL_K->GetNbinsX());
  }
  
  return wl;
}

//Apply snells law as if we were refracting from Si->SiO2->LXe
//If we're greater than the critical angle, the photons are "reflected" back into the SiPM
void ExternalCrossTalkProcess::ApplySnellsLaw(double wl, ROOT::Math::XYZVector &dir, ROOT::Math::XYZVector &norm)
{

  //This may be slow but we can update later
  G4double n_lxe = gLXeRefractiveWavelength->Eval(wl);
  G4double n_si = gSiRefractiveWavelength->Eval(wl);
  G4double n_sio2 = gSiO2RefractiveWavelength->Eval(wl);
  
  //G4cout<<" WL = " << wl << " n_lxe = "<< n_lxe << " n_si = "<< n_si << " n_sio2 = "<< n_sio2 <<G4endl;

  G4double ndoti = dir.Dot(norm);//Incident dot normal
  G4double theta = std::acos(ndoti);//Angle wrt incidence

  //Calculated refracted angle wrt surface normal using snells law
  //Std methods catch illegal angles and return NaN
  theta = std::asin(n_si/n_sio2*std::sin(theta)); //Refract from Si to SiO2

  //What do if it is nan? Kill instead of internal reflection to avoid double-counting direct CT
  if(std::isnan(theta)){
    //dir = dir - (2.*ndoti)*norm;//Perfectly specular reflection back into SiPM
    dir = ROOT::Math::XYZVector(0.0,0.0,-100.0);//Flag to kill photon by setting dir to -100.0
  } else { 
    G4double mu = n_si/n_sio2;
    //Snells law in 3D  
    ROOT::Math::XYZVector rVR = std::sqrt(1-mu*mu*(1-ndoti*ndoti))*norm + mu*(dir-ndoti*norm);
    mu = n_sio2/n_lxe;//Apply next interface
    ndoti = rVR.Dot(norm);
    theta = std::acos(ndoti);//Get new angle wrt to surface normal
    theta = std::asin(n_sio2/n_lxe*std::sin(theta));//Snells again
    if(std::isnan(theta)){
      //dir = dir - (2.*ndoti)*norm;//Perfectly specular reflection back into SiPM
       dir = ROOT::Math::XYZVector(0.0,0.0,-100.0);//Flag to kill photon by setting dir to -100.0
    }else{
      rVR = std::sqrt(1-mu*mu*(1-ndoti*ndoti))*norm + mu*(rVR-ndoti*norm);//Snells law applied in 3D
      dir = rVR;
      dir = dir.Unit();//Normalize to unit
    }
  }

} 
///Update May 3rd, the refraction code here is bugged, not used anymore, left for historical purposes
//Apply simulation input as if we were refracting from Si->SiO2->LXe
//Uses Kurtis' ANSYS Lumerical far field EM simulation for cross-talk propogation from within the SiPM to LXe
void ExternalCrossTalkProcess::ApplySimInput(double wl, ROOT::Math::XYZVector &dir, ROOT::Math::XYZVector &norm)
{
  double wl_um = wl/1000.0;//Data in hTransmission and hConversion is in um
  G4double ndoti = dir.Dot(norm);//Incident dot normal
  G4double theta = std::acos(ndoti);//Angle from surface normal to incidence ray
  theta = 180.0*theta/TMath::Pi(); //Convert to degrees
  double prob_trans = hTransmission->GetBinContent(hTransmission->FindBin(wl_um,theta));
  double p_rng = gRandom->Uniform();//Random number for probability throw
  if(p_rng < prob_trans){
    //Transform Dir into new direction
    //Get data from hConversion for new theta
    double theta_refr = hConversion->GetBinContent(hConversion->FindBin(wl_um,theta));// "Refract" into LXe from Si
    double diff_theta = theta_refr-theta;
    //G4cout << "Theta before = " << theta << " theta after = "<< theta_refr << " difference = "<< diff_theta << G4endl;
    theta_refr = TMath::Pi()*theta_refr/180.0;//Convert back to radians
    //Project down to "refract" the photon path
    double projZ = sqrt(dir.Mag2())*sin(diff_theta);
    dir = dir - ROOT::Math::XYZVector(0,0,projZ);
    dir = dir.Unit();
  } else {
    dir = ROOT::Math::XYZVector(0.0,0.0,-100.0);//Flag to kill photon by setting dir to -100.0
  }
}

//Apply simulation input as if we were refracting from Si->SiO2->LXe
//Uses Kurtis' ANSYS Lumerical far field EM simulation for cross-talk propogation from within the SiPM to LXe
//Manual refraction using Snells law twice instead of conversion table
void ExternalCrossTalkProcess::ApplySimInputSnells(double wl, ROOT::Math::XYZVector &dir, ROOT::Math::XYZVector &norm)
{
//Get Refractive indices
//This may be slow but we can update later
  G4double n_lxe = gLXeRefractiveWavelength->Eval(wl);
  G4double n_si = gSiRefractiveWavelength->Eval(wl);
  G4double n_sio2 = gSiO2RefractiveWavelength->Eval(wl);
  // G4cout << Form("%f,%f,%f \n",n_si,n_sio2,n_lxe);

  G4double ndoti = dir.Dot(norm);//Incident dot normal
  G4double theta = std::acos(ndoti);//Angle from surface normal to incidence ray
  // G4cout << "Internal Theta "<<theta/TMath::Pi()*180.0 <<G4endl;
  // G4cout << "Internal ThetaDir "<<dir.theta()/TMath::Pi()*180.0 <<G4endl;
  G4double wl_um = wl/1000.0;
  G4double prob_trans = hTransmission->GetBinContent(hTransmission->FindBin(wl_um,theta*180.0/TMath::Pi()));
  G4double p_rng = gRandom->Uniform();//Random number for probability throw

  if(p_rng < prob_trans){
    //Transform Dir into new direction
    //Get data from hConversion for new theta
    // cout << "Theta internal "<< theta/TMath::Pi()*180.0<<endl;
    // double theta_r = hConversion->GetBinContent(hConversion->FindBin(wl_um,theta/TMath::Pi()*180.0));// "Refract" into LXe from Si
    // G4cout << "Theta Refracted = "<< theta_r <<G4endl;

    //Snells law
    //Calculated refracted angle wrt surface normal
    theta = std::asin(n_si/n_sio2*std::sin(theta)); 
    // G4cout << Form("Dir before = [%f,%f,%f]\n",dir.X(),dir.Y(),dir.Z());
    // G4cout << "Theta from vector " <<dir.theta()/TMath::Pi()*180.0 <<G4endl;
    double mu = n_si/n_sio2;
    //Snells law in 3D  
    ROOT::Math::XYZVector rVR = std::sqrt(1-mu*mu*(1-ndoti*ndoti))*norm + mu*(dir-ndoti*norm);
    mu = n_sio2/n_lxe;//Apply next interface
    ndoti = rVR.Dot(norm);
    //Refract again into lxe
    theta = std::acos(ndoti);
    // G4cout << "Theta after SiO2 = "<< theta*180.0/TMath::Pi() <<G4endl;
    theta = std::asin(mu*std::sin(theta));
    // G4cout << "Theta after both = "<< theta*180.0/TMath::Pi() <<G4endl;
    rVR = std::sqrt(1-mu*mu*(1-ndoti*ndoti))*norm + mu*(rVR-ndoti*norm);//Snells law applied in 3D
    dir = rVR;
    dir = dir.Unit();//Normalize to unit
    ndoti = rVR.Dot(norm);
    //Refract again into lxe
    theta = std::acos(ndoti);
    // G4cout << Form("Dir after = [%f,%f,%f]\n",dir.X(),dir.Y(),dir.Z());
    // G4cout << "Theta from Vector = "<< theta*180.0/TMath::Pi() <<G4endl;
   } else {
    dir = ROOT::Math::XYZVector(0.0,0.0,-100.0);//Flag to kill photon by setting dir to -100.0
  }
}

//Read ROOT Files from data dir
void  ExternalCrossTalkProcess::ReadROOTFiles()
{ 
  G4cout << "ExternalCrossTalkProcess: Reading ROOT input files"<< G4endl;
  inputFile = new TFile(Form("%s/%s",std::getenv("LOLXSTLDIR"),"lolx_sim_data.root"),"READ");
  gLXeRefractiveWavelength = new TGraph((TH1D*) inputFile->Get("rindex_lxe"));
  gSiRefractiveWavelength = new TGraph((TH1D*) inputFile->Get("rindex_si"));
  gSiO2RefractiveWavelength = new TGraph((TH1D*) inputFile->Get("rindex_sio2"));
  hEXTWavelength = (TH1D*) inputFile->Get("ext_wavelength");

  inputFileKurtis = new TFile(Form("%s/%s",std::getenv("LOLXSTLDIR"),"HPK_140nm_Production_And_Transmission_highd.root"),"READ");
  hEXTWL_K = (TH1D*) inputFileKurtis->Get("ProdCurve");
  hTransmission = (TH2D*) inputFileKurtis->Get("SiO2Trans");//Probability to transmit from Si->SiO2
  hConversion = (TH2D*) inputFileKurtis->Get("AngleConv");// Angle in LXe (Z-axis) for 
}
