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
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  DefineMaterials();

  return DefineVolumes();

}
void DetectorConstruction::DefineMaterials()
{
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Fe");
  nistManager->FindOrBuildMaterial("G4_lAr");

  G4double z;
  G4double a;
  G4double density;
  new G4Material("Galactic", z = 1., a = 1.01 * g / mole, density = universe_mean_density, kStateGas, 2.73 * kelvin, 3.e-18 * pascal);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  //Liquid Argon Volume:
  //https://cds.cern.ch/record/2884152/files/2312.09104.pdf, Table 2.1

  //Primary steel membrane: d = 1.2mm
  //https://arxiv.org/pdf/2312.09104, Table 1.2
  G4double height = 14.0 * m;
  G4double width = 15.1 * m;
  G4double length = 62.0 * m;
  G4double primaryIronThickness = 1.2 * mm; //1.2 * mm
  G4double PUThickness = 80.0 * cm;
  G4double outerIronThickness = 54.63* mm; //54.63* mm see DUNE_FD_analysis/helper.py for calculation
  G4double cryostatThickness = primaryIronThickness + PUThickness + outerIronThickness;

  auto xPrimaryIron = length + 2 * primaryIronThickness;
  auto yPrimaryIron = width + 2 * primaryIronThickness;
  auto zPrimaryIron = height + 2 * primaryIronThickness;

  auto xPU = xPrimaryIron + 2* PUThickness;
  auto yPU = yPrimaryIron + 2* PUThickness;
  auto zPU = zPrimaryIron + 2* PUThickness;

  auto xOuterIron = xPU + 2* outerIronThickness;
  auto yOuterIron = yPU + 2* outerIronThickness;
  auto zOuterIron = zPU + 2* outerIronThickness;

  auto xTotal = length + 2 * cryostatThickness;
  auto yTotal = width + 2 * cryostatThickness;
  auto zTotal = height + 2 * cryostatThickness;

  auto worldSizeX = 1.2 * xTotal;
  auto worldSizeY = 1.2 * yTotal;
  auto worldSizeZ = 1.2 * zTotal;

  auto* nist = G4NistManager::Instance();

  auto* PUMaterial = new G4Material(
      "PUFoam",
      0.090*g/cm3,   // 90 kg/m^3
      4               
  );

  // Mass fractions from https://pubs.rsc.org/an/article/149/5/1609/795207/Toward-chemical-recycling-of-PU-foams-study-of-the
  PUMaterial->AddElement(nist->FindOrBuildElement("H"), 0.0660);
  PUMaterial->AddElement(nist->FindOrBuildElement("C"), 0.6780);
  PUMaterial->AddElement(nist->FindOrBuildElement("N"), 0.0763);
  PUMaterial->AddElement(nist->FindOrBuildElement("O"), 0.1797);

  
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto ironMaterial = G4Material::GetMaterial("G4_Fe");
  auto argonMaterial = G4Material::GetMaterial("G4_lAr");

  auto worldS = new G4Box("World",  // its name
                          worldSizeX / 2, worldSizeY / 2, worldSizeZ / 2);  // its size

  auto worldLV = new G4LogicalVolume(worldS,  // its solid
                                     defaultMaterial,  // its material
                                     "World");  // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
                                   G4ThreeVector(),  // at (0,0,0)
                                   worldLV,  // its logical volume
                                   "World",  // its name
                                   nullptr,  // its mother  volume
                                   false,  // no boolean operation
                                   0,  // copy number
                                   fCheckOverlaps);  // checking overlaps

  auto outerIronS = new G4Box("OuterIron", xOuterIron / 2, yOuterIron / 2, zOuterIron/ 2);
  auto outerIronLV = new G4LogicalVolume(outerIronS, ironMaterial, "OuterIron");
  new G4PVPlacement(nullptr, G4ThreeVector(), outerIronLV, "OuterIron", worldLV, false, 0, fCheckOverlaps);
  
  auto PUS = new G4Box("PU", xPU / 2, yPU / 2, zPU/ 2);
  auto PULV = new G4LogicalVolume(PUS, PUMaterial, "PU");
  new G4PVPlacement(nullptr, G4ThreeVector(), PULV, "PU", outerIronLV, false, 0, fCheckOverlaps);


  auto primaryIronS = new G4Box("PrimaryIron", xPrimaryIron / 2, yPrimaryIron / 2, zPrimaryIron / 2);
  auto primaryIronLV = new G4LogicalVolume(primaryIronS, ironMaterial, "PrimaryIron");
  new G4PVPlacement(nullptr, G4ThreeVector(), primaryIronLV, "PrimaryIron", PULV, false, 0, fCheckOverlaps);

  auto argonS = new G4Box("lAr", length / 2, width / 2, height / 2);
  auto argonLV = new G4LogicalVolume(argonS, argonMaterial, "lAr");

  fScoringVolume = argonLV;
  //auto argonVis = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
  //argonVis->SetForceWireframe(true);
  //argonLV->SetVisAttributes(argonVis);
  
  new G4PVPlacement(nullptr, G4ThreeVector(), argonLV, "lAr", primaryIronLV, false, 0, fCheckOverlaps);


  return worldPV;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
