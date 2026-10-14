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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B1::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <cmath>
#include <numbers>
#include <iostream>
#include "Randomize.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fMessenger = new G4GenericMessenger(this,"/FDSteel/","FDSteel generator controls");

  fMessenger->DeclarePropertyWithUnit("electronEnergy", "MeV", fElectronEnergy, "Set electron kinetic energy");

  //fMessenger->DeclarePropertyWithUnit("nuEnergy","MeV",nuEnergy,"Set true neutrino energy");
  //fMessenger->DeclareProperty("channel", channel,"Set interaction channel");
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);
}
G4ThreeVector PrimaryGeneratorAction::SampleFePosition()
{
  G4double L = 62.0 * m;
  G4double W = 15.1 * m;
  G4double H = 14.0 * m;
  G4double d = 62.7 * mm;

  G4double Lout = L + 2.0 * d;
  G4double Wout = W + 2.0 * d;
  G4double Hout = H + 2.0 * d;

  G4double VFrontBack = Lout * Hout * d;
  G4double VLeftRight = d * W * Hout;
  G4double VTopBottom = d * L * W;

  G4double Vtotal = 2 * VFrontBack + 2 * VLeftRight + 2 * VTopBottom;


  G4double xpos;
  G4double ypos;
  G4double zpos; 

  G4double random = G4UniformRand();
  G4double FB_threshold = 2 * VFrontBack / Vtotal;
  G4double LR_threshold = 2 * VLeftRight / Vtotal;
  //G4double TB_threshold = 2 * VTopBottom / Vtotal;



  G4double which_plate = G4UniformRand();
  if (random < FB_threshold)
  {
    
    if (which_plate < 0.5)
    {
      //FRONT
      xpos = (G4UniformRand() - 0.5) * Lout;
      ypos = W / 2 + (G4UniformRand() * d);
      zpos = (G4UniformRand() - 0.5) * Hout;
      
    }
    else
    {
      //BACK
      xpos = (G4UniformRand() - 0.5) * Lout;
      ypos = - W / 2 - (G4UniformRand() * d) ;
      zpos = (G4UniformRand() - 0.5) * Hout;
      
    }
  }
  
  else if ((random >= FB_threshold) && (random < FB_threshold + LR_threshold))
  {
    if (which_plate < 0.5)
    {
      //LEFT
      xpos = L / 2 + (G4UniformRand() * d);
      ypos = (G4UniformRand() - 0.5) * W;
      zpos = (G4UniformRand() - 0.5) * Hout;
    }
    else
    {
      //RIGHT
      xpos = - L / 2 - (G4UniformRand() * d) ;
      ypos = (G4UniformRand() - 0.5) * W;
      zpos = (G4UniformRand() - 0.5) * Hout;
    }
  }
  else
  {
    if (which_plate < 0.5)
    {
      //TOP
      xpos = (G4UniformRand() - 0.5) * L;
      ypos = (G4UniformRand() -0.5) * W;
      zpos = H / 2 + (G4UniformRand() * d);
    }
    else
    {
      //BOTTOM
      xpos = (G4UniformRand() - 0.5) * L;
      ypos = (G4UniformRand() -0.5) * W;
      zpos = - H / 2 - (G4UniformRand() * d);
    }
  }
  return G4ThreeVector(xpos, ypos, zpos);

}

G4ThreeVector PrimaryGeneratorAction::SampleIsotropicDirection()
{
  G4double u = G4UniformRand();
  G4double v = G4UniformRand();
  G4double phi = CLHEP::twopi * u;
  G4double cos_theta = 2*v - 1;
  G4double sin_theta = std::sqrt(1 - cos_theta * cos_theta);

  G4double x = sin_theta * std::cos(phi);
  G4double y = sin_theta * std::sin(phi);
  G4double z = cos_theta;
  return G4ThreeVector(x, y, z);
}

void PrimaryGeneratorAction::GenerateElectron(G4Event* event)
{
  auto particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  fParticleGun->SetParticleDefinition(particle);
  //G4double energy = 15 * MeV;

  fParticleGun->SetParticleEnergy(fElectronEnergy);
  fParticleGun->SetParticlePosition(SampleFePosition());
  fParticleGun->SetParticleMomentumDirection(SampleIsotropicDirection());

  fParticleGun->GeneratePrimaryVertex(event);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  GenerateElectron(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
