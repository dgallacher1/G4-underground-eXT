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
/// \file optical/LXe/src/LXeMainVolume.cc
/// \brief Implementation of the LXeMainVolume class
//
//
#include "VUV4Volume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double VUV4Volume::Package_sizeXY = 15*mm;  // size of the  MPPC package initially made from Xenon.
G4double VUV4Volume::Package_sizeZ = 2.5*mm;  // everything is housed within this volume
G4double VUV4Volume::Package_border = 1.1*mm;


VUV4Volume::VUV4Volume()
  //Pass info to the G4PVPlacement constructor
  :G4LogicalVolume(new G4Box("Package", Package_sizeXY/2, Package_sizeXY/2, Package_sizeZ/2),DMaterials::GetDetectorMedium(), "Package")
{
    G4cout <<"Building VUV4Volume"<<G4endl;
    bool checkOverlaps = true;

    // create a logical volume, where will put everything

    G4double SiPM_sizeXY = 5.9*mm; // size of individual MPPC, which is 5.85 x 5.95 in actuallity (approximated by a square of 5.9x5.9 here)
    G4double SiPM_sizeZ = 0.05*mm; // thickness of MPPC
    // position of MPPC arrays inside the package (Packagelogic)

    SiPM_planeG4 = G4ThreeVector(0.0,0.0,SiPM_sizeZ/2.);
    SiPM_plane = TVector3(0.0,0.0,SiPM_sizeZ/2.);
    //MPPC
    G4Box* SolidSiPM = new G4Box("SiPMsold", SiPM_sizeXY/2., SiPM_sizeXY/2., SiPM_sizeZ/2.);
    SiPMlogic = new G4LogicalVolume(SolidSiPM, DMaterials::Get_silicon_mat(), "SiPMlog"); // MPPC, made out of silicon

    // To form a ceramic holding structure need to do a boolean subtraction of volumes
    // first creating two boxes, one the size of the Packagelogic, the other is half depth (z) and is smaller by the Package_border

    G4double delta = 0.0*mm;

    G4Box* box1 = new G4Box("Box1",Package_sizeXY/2.,Package_sizeXY/2.,Package_sizeZ/2.); // size of the Packagelogic
    G4Box* box2 = new G4Box("Box2",(Package_sizeXY-Package_border)/2.,(Package_sizeXY-Package_border)/2.,Package_sizeZ/4.+delta); // half depth
                                                                                               // and smaller in size by Package_border
    G4ThreeVector  translation_ceramic(0.0,0.0,(Package_sizeZ/4.)+delta/2.); // this vectro is needed for volume subtraction
    G4SubtractionSolid* b1minusb2 = new G4SubtractionSolid("box1-box2",box1,box2,0,translation_ceramic); // creating a solid for ceramic holding structure
    //by subtracting one G4Box from another G4Box.
    //Ceramiclogic = new G4LogicalVolume(b1minusb2, DMaterials::Get_silicon_mat(), "ceramic"); // logical volume for the ceramic holding structure
    //Ceramiclogic = new G4LogicalVolume(box1, DMaterials::Get_xenon_mat(), "ceramic");
    Ceramiclogic = new G4LogicalVolume(b1minusb2, DMaterials::Get_fCeramic(), "ceramic");
    new G4PVPlacement(0, G4ThreeVector(), Ceramiclogic, "Ceramic", this,0, checkOverlaps);

    // creating a quarts window
    G4double Window_sizeZ = 0.5*mm; // window thickness  nominally was .5mm
    G4double Window_gap_scaler = 1.; // normall 1.5 how much larger the distance from the edge of the Packagelogic will be compared to the Package_border
    G4double Window_sizeXY = Package_sizeXY-Window_gap_scaler*Package_border;

    G4Box* window_solid = new G4Box("QuartzWindow",Window_sizeXY/2.,Window_sizeXY/2.,Window_sizeZ/2.);
    Windowlogic = new G4LogicalVolume(window_solid, DMaterials::Get_quartz_mat(), "QuartzWindow");
    //Add a small buffer volume between window and SiO2 to capture photons
    G4double buffer_halfheight = 0.05*mm;
    G4Box* buffer_solid = new G4Box("LXeBufferBox",Window_sizeXY/2.,Window_sizeXY/2.,buffer_halfheight);
    BufferLogic = new G4LogicalVolume(buffer_solid, DMaterials::GetDetectorMedium(), "LXeBuffer");

    SiliconPMs.push_back(new G4PVPlacement(0, SiPMProperties::GetChipOffsetG4(0)+SiPM_planeG4, SiPMlogic, "SiPM0", this,false,0, checkOverlaps)); // individual MPPCs
    SiliconPMs.push_back(new G4PVPlacement(0, SiPMProperties::GetChipOffsetG4(1)+SiPM_planeG4, SiPMlogic, "SiPM1", this,false,1, checkOverlaps));
    SiliconPMs.push_back(new G4PVPlacement(0, SiPMProperties::GetChipOffsetG4(2)+SiPM_planeG4, SiPMlogic, "SiPM2", this,false,2, checkOverlaps));
    SiliconPMs.push_back(new G4PVPlacement(0, SiPMProperties::GetChipOffsetG4(3)+SiPM_planeG4, SiPMlogic, "SiPM3", this,false,3, checkOverlaps));

    G4double SiO2_thickness = 110.0*nm;
    G4double window_offset = 0.5*mm;
    //Shift window up to make room
    WindowFaceZ = Window_sizeZ/2.+SiPM_sizeZ+window_offset;
    new G4PVPlacement(0, G4ThreeVector(0.0,0.0,WindowFaceZ), Windowlogic, "QuartzWindow", this,0, checkOverlaps);
    //Buffer volume between window and SiO2 layers
    LXeBuffer = new G4PVPlacement(0, G4ThreeVector(0.0,0.0,SiPM_sizeZ+SiO2_thickness+buffer_halfheight), BufferLogic, "LXeBuffer", this,0, checkOverlaps);
    //SiO2 Layers
    G4ThreeVector SiO2Plane(0.0,0.0,SiPM_sizeZ+SiO2_thickness/2.0);
    G4Box* SolidSiO2 = new G4Box("BoxSiO2",SiPM_sizeXY/2.,SiPM_sizeXY/2.,SiO2_thickness/2.);
    SiO2Logic = new G4LogicalVolume(SolidSiO2, DMaterials::Get_quartz_mat(), "SiO2log"); //made out of sio2
    SiO2Layers.push_back(new G4PVPlacement(0, SiPMProperties::GetChipOffsetG4(0)+SiO2Plane, SiO2Logic, "SiO20", this,false,0, checkOverlaps)); // individual SiO2 layers
    SiO2Layers.push_back(new G4PVPlacement(0, SiPMProperties::GetChipOffsetG4(1)+SiO2Plane, SiO2Logic, "SiO21", this,false,1, checkOverlaps));
    SiO2Layers.push_back(new G4PVPlacement(0, SiPMProperties::GetChipOffsetG4(2)+SiO2Plane, SiO2Logic, "SiO22", this,false,2, checkOverlaps));
    SiO2Layers.push_back(new G4PVPlacement(0, SiPMProperties::GetChipOffsetG4(3)+SiO2Plane, SiO2Logic, "SiO23", this,false,3, checkOverlaps));
    SurfaceProperties();
}

G4ThreeVector VUV4Volume::GetChipPosition(G4int num,G4RotationMatrix* rot,G4ThreeVector trans){

    TVector3 SiPM_pos = SiPMProperties::GetChipOffset(num)+SiPM_plane;
    G4double x = SiPM_pos.x()*rot->xx()+SiPM_pos.y()*rot->yx()+SiPM_pos.z()*rot->zx()+trans.x();
    G4double y = SiPM_pos.x()*rot->xy()+SiPM_pos.y()*rot->yy()+SiPM_pos.z()*rot->zy()+trans.y();
    G4double z = SiPM_pos.x()*rot->xz()+SiPM_pos.y()*rot->yz()+SiPM_pos.z()*rot->zz()+trans.z();
    return G4ThreeVector(x,y,z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VUV4Volume::SurfaceProperties(){

    //Flag for systematic studies
    int efficiency_flag = 0; // 0 == nominal, 1 == 10% rel. increase above 600nm, 2 == 10% rel. decrease above 600nm

    const G4int NUMENTRIES1 = 850;  //from 150 to 1000
    G4double energy[NUMENTRIES1];
    G4double reflectivity[NUMENTRIES1];
    G4double efficiency[NUMENTRIES1];
    G4double reflectivity_ceramic[NUMENTRIES1];
    G4double efficiency_ceramic[NUMENTRIES1];
    G4double reflectivity_window[NUMENTRIES1];

    for(int iE=0; iE < NUMENTRIES1; iE++){
      G4double WL=1000.-(iE); // starting from 1000nm WL, ending with 150nm
      energy[iE]=1240./WL; // tabulate energy
      reflectivity[iE] = 0.0;//Set reflectivity to 0 for all WLs, no longer used -DG
      //Add reflectivity for photons above 600 nm (EXT)
      //efficiency[iE] = 1.0;//100% efficiency to check reflectivity
      efficiency[iE] =  LOLXReadData::GetSiPM_Efficiency(energy[iE]); // Temporarily set PDE to unity for all WLs -DG
      //G4cout << Form("SiPM Eff (eV,E) = (%f,%f)\n",energy[iE],efficiency[iE]);
      if(WL > 600) {
        if(efficiency_flag == 1) efficiency[iE] = efficiency[iE]*1.1;
        else if(efficiency_flag == 2) efficiency[iE] = efficiency[iE]*0.9;
      }

      reflectivity_ceramic[iE] = 0.25;
      efficiency_ceramic[iE] = 0.0;
       // printf("Wl = %f energy = %f eff = %f\n",WL,energy[iE],efficiency[iE]);
      energy[iE] *= eV;
      reflectivity_window[iE] = 0.05;//No longer used - DG
    }

    //Add SiPM Optics- DG
    //Get Si Data
    std::vector<double> si_real_vec;
    std::vector<double> si_wl_vec;
    std::vector<double> si_imag_vec;

    GetSiData(si_wl_vec,si_real_vec,si_imag_vec);
    G4double* silicon_real_index = (G4double*) si_real_vec.data();//Static cast from float* to double*
    G4double* silicon_imaginary_index = (G4double*) si_imag_vec.data();
    G4double* silicon_wl = (G4double*) si_wl_vec.data();
    G4int data_length = si_wl_vec.size();
    
    //From Austin, perfect transmission so G4 does reflectivity calculations properly -- Causes things to crash, maybe due to different G4 versions from Austin to our sim
    G4double t = 0.0;
    G4double photonE_noTransmit[2] = {1.0*eV, 12.4*eV};
    G4double noTransmit[2] = {t, t};
    G4double noRef[2] = {1.0-t, 1.0-t};
       
    LOLXOpticalSurface* OpSiPMSurface = new LOLXOpticalSurface("OpSiPMSurface");
    OpSiPMSurface->SetModel(lunified); //unified
    OpSiPMSurface->SetType(dielectric_metal); //_metal
    OpSiPMSurface->SetFinish(lpolished);
    OpSiPMSurface->SetSigmaAlpha(0.0);

    G4MaterialPropertiesTable* SiPMTable = new G4MaterialPropertiesTable();
    SiPMTable->AddProperty("REALRINDEX",silicon_wl,silicon_real_index,data_length);
    SiPMTable->AddProperty("TRANSMITTANCE",photonE_noTransmit,noTransmit,2);
    SiPMTable->AddProperty("IMAGINARYRINDEX",silicon_wl,silicon_imaginary_index,data_length);
    // SiPMTable->AddProperty("REFLECTIVITY",energy,reflectivity,NUMENTRIES1); // Disable reflectivity - If not present it will calculate from Fresnel using Real and Imag index (See LOLXMaterials.hh)
    SiPMTable->AddProperty("EFFICIENCY",energy,efficiency,NUMENTRIES1);
    OpSiPMSurface->SetMaterialPropertiesTable(SiPMTable);
    new G4LogicalSkinSurface("OpSiPMSurface",SiPMlogic,OpSiPMSurface);

    //SiO2/Quartz surface 
    LOLXOpticalSurface* WindowSurface = new LOLXOpticalSurface("WindowSurface");
    WindowSurface->SetType(dielectric_dielectric);
    WindowSurface->SetModel(lunified);
    WindowSurface->SetFinish(lpolished);
    WindowSurface->SetSigmaAlpha(0.);

    //Get SiO2 Data
    std::vector<double> sio2_real_vec;
    std::vector<double> sio2_wl_vec;
    std::vector<double> sio2_imag_vec;

    GetSiO2Data(sio2_wl_vec,sio2_real_vec,sio2_imag_vec);
    G4double* sio2_real_index = (G4double*) sio2_real_vec.data();
    G4double* sio2_imaginary_index = (G4double*) sio2_imag_vec.data();
    G4double* sio2_wl = (G4double*) sio2_wl_vec.data();
    data_length = sio2_wl_vec.size();

    G4MaterialPropertiesTable* WindowTable = new G4MaterialPropertiesTable();
    WindowTable->AddProperty("RINDEX", sio2_wl, sio2_real_index,data_length);
    WindowSurface->SetMaterialPropertiesTable(WindowTable);
    new G4LogicalSkinSurface("WindowSurface",Windowlogic,WindowSurface);

    //Add SiO2 Boundary Surfaces
    LOLXOpticalSurface* SiO2Surface = new LOLXOpticalSurface("SiO2Surface");
    SiO2Surface->SetModel(lunified); //unified
    SiO2Surface->SetType(dielectric_dielectric); 
    SiO2Surface->SetFinish(lpolished);
    SiO2Surface->SetMaterialPropertiesTable(WindowTable);//Same as window
    //Change skin on SiO2 to border->In for reflectivity and absorption out (In MainVolume_Cylindrical)
    // new G4LogicalBorderSurface("inSiO2Surface0",LXeBuffer,SiO2Layers[0],WindowSurface);
    // new G4LogicalBorderSurface("inSiO2Surface1",LXeBuffer,SiO2Layers[1],WindowSurface);
    // new G4LogicalBorderSurface("inSiO2Surface2",LXeBuffer,SiO2Layers[2],WindowSurface);
    // new G4LogicalBorderSurface("inSiO2Surface3",LXeBuffer,SiO2Layers[3],WindowSurface);

    LOLXOpticalSurface* OpCeramicSurface = new LOLXOpticalSurface("OpCeramicSurface");
    OpCeramicSurface->SetModel(lglisur);
    OpCeramicSurface->SetType(dielectric_metal);
    OpCeramicSurface->SetFinish(lground);

    for(double &r: reflectivity_ceramic) r=0.0; // 75% Reflective ceramic MPPC holding structure
    for(double &r: efficiency_ceramic) r=0; // zero efficiency
    G4MaterialPropertiesTable* CeramicTable = new G4MaterialPropertiesTable();
    CeramicTable->AddProperty("REFLECTIVITY",energy,reflectivity_ceramic,NUMENTRIES1);
    CeramicTable->AddProperty("EFFICIENCY",energy,efficiency_ceramic,NUMENTRIES1);
    OpCeramicSurface->SetMaterialPropertiesTable(CeramicTable);
    new G4LogicalSkinSurface("OpCeramicSurface",Ceramiclogic,OpCeramicSurface); //This surface is for the Ceramic, hence for Ceramiclogic

}

//Read in silicon data from Austin - DG
void VUV4Volume::GetSiData(std::vector<double>&wl, std::vector<double>&real_index, std::vector<double>&k_index){
  const std::string filename_real = Form("%s/%s",std::getenv("LOLXSTLDIR"),"/r_index/silicon_merged_n_ev.txt");
  std::ifstream file(filename_real);
  //Read real index values
  // G4cout << "Reading Si R-Index data.."<<G4endl;

  if (file.is_open()) {
      std::string line;
      // Skip the first line
      std::getline(file, line);
      // G4cout << "Real data [(n,eV),(),...]"<<G4endl;
      // G4cout << "[";
      while (std::getline(file, line)) {
          std::istringstream iss(line);//input string stream
          float val1, val2;
          if (iss >> val1 >> val2) {
              wl.push_back(val1*eV);
              real_index.push_back(val2);
              // G4cout << Form("(%f,%f),",val1*eV,val2);
          }
      }
      // G4cout << "]"<<G4endl;
      file.close();
  } else {
      std::cerr << "Unable to open file: " << filename_real << std::endl;
  }

  const std::string filename_imag =  Form("%s/%s",std::getenv("LOLXSTLDIR"),"/r_index/silicon_merged_k_ev.txt");
  std::ifstream file_imag(filename_imag);
  //Read coefficient of extinction values
  if (file_imag.is_open()) {
      std::string line;
      // Skip the first line
      std::getline(file_imag, line);
      while (std::getline(file_imag, line)) {
          std::istringstream iss(line);//input string stream
          float val1, val2;
          if (iss >> val1 >> val2) {
              k_index.push_back(val2);
          }
      }
      file_imag.close();
  } else {
      std::cerr << "Unable to open file: " << filename_imag << std::endl;
  }

}
//Read in silicon data from Austin - DG
void VUV4Volume::GetSiO2Data(std::vector<double>&wl, std::vector<double>&real_index, std::vector<double>&k_index){
  const std::string filename_real = Form("%s/%s",std::getenv("LOLXSTLDIR"),"/r_index/sio2_n_eV.txt");
  std::ifstream file(filename_real);
  //Read real index values
//   G4cout << "Reading SiO2 R-Index data.."<<G4endl;

  if (file.is_open()) {
      std::string line;
      // Skip the first line
      std::getline(file, line);
    //   G4cout << "Real data [(n,eV),(),...]"<<G4endl;
    //   G4cout << "[";
      while (std::getline(file, line)) {
          std::istringstream iss(line);//input string stream
          float val1, val2;
          if (iss >> val1 >> val2) {
              wl.push_back(val1*eV);
              real_index.push_back(val2);
            //   G4cout << Form("(%f,%f),",val1*eV,val2);
          }
      }
    //   G4cout << "]"<<G4endl;
      file.close();
  } else {
      std::cerr << "Unable to open file: " << filename_real << std::endl;
  } 


  const std::string filename_imag = Form("%s/%s",std::getenv("LOLXSTLDIR"),"/r_index/sio2_k_eV.txt");
  std::ifstream file_imag(filename_imag);
  //Read coefficient of extinction values
  if (file_imag.is_open()) {
      std::string line;
      // Skip the first line
      std::getline(file_imag, line);
      while (std::getline(file_imag, line)) {
          std::istringstream iss(line);//input string stream
          float val1, val2;
          if (iss >> val1 >> val2) {
              k_index.push_back(val2);
          }
      }
      file_imag.close();
  } else {
      std::cerr << "Unable to open file: " << filename_imag << std::endl;
  }

}
