/// @file
/// @copyright Copyright (C) 2016-2026 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

// Standard library includes
#include <iostream>
#include <array>
#include <cmath>

// Geant4 includes
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
//#include "globals.hh"

// MARLEY includes
#include "marley/Error.hh"
#include "marley/JSONConfig.hh"
#include "marley/hepmc3_utils.hh"

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

// marg4 includes
#include "MarleyPrimaryGeneratorAction.hh"


namespace B1
{


MarleyPrimaryGeneratorAction::MarleyPrimaryGeneratorAction(
  const std::string& config_file_name) : G4VUserPrimaryGeneratorAction()
{
  // Create a new marley::Generator object using the settings from the
  // configuration file.
  
fEnu = 40.0 * MeV;
fMessenger = std::make_unique<G4GenericMessenger>(this,"/marley/","MARLEY primary generator control");

fMessenger->DeclarePropertyWithUnit("nuEnergy","MeV",fEnu,"Set the incoming neutrino energy.");
fMessenger->DeclareProperty("ironVolume",ironVolume,"Set the generation iron volume.");
fMessenger->DeclareProperty("nuFlavor",nuFlavor,"Set the incoming neutrino flavor.");
  marley::JSONConfig config( config_file_name );

  marley_generator_= config.create_generator();
}

void MarleyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Create a new primary vertex at the spacetime origin.
  G4ThreeVector position;

  if (ironVolume == "PrimaryIron")
  {
    position = SamplePrimaryIronPosition();
  }
  else if (ironVolume == "OuterIron")
  {
    position = SampleOuterIronPosition();
  }

  else
  {
      G4Exception(
          "MarleyPrimaryGeneratorAction::GeneratePrimaries",
          "InvalidIronVolume",
          FatalException,
          "Unknown iron generation volume!"
      );
      return;
  }

  G4PrimaryVertex* vertex = new G4PrimaryVertex(position.x(), position.y(), position.z(), 0.); // x,y,z,t0
 
  // Generate a new MARLEY event using the owned marley::Generator object
   // MeV
  std::array<double, 3> nuDir = SampleIsotropicDirection();
  


  auto ev = marley_generator_.create_event(nuFlavor,fEnu / MeV,1000260560,nuDir);

  // Account for possibly different systems of units for the particle
  // 4-momenta by using this conversion factor
  auto mom4_conv_factor = MeV;

  // Query the event to determine what it uses for energy/momentum units
  auto ev_energy_unit = ev->momentum_unit();
  if ( ev_energy_unit == HepMC3::Units::GEV ) {
    mom4_conv_factor = GeV;
  }
  else if ( ev_energy_unit != HepMC3::Units::MEV ) {
    throw marley::Error( "Unrecognized momentum unit encountered in"
      " MarleyPrimaryGeneratorAction::GeneratePrimaries()" );
  }

  // Collect the final-state particles from the HepMC3 event
  auto finals = marley_hepmc3::get_particles_with_status(
    marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, *ev);

  // Loop over each of the final particles in the MARLEY event
  for ( const auto& fp : finals ) {

    // Access the 4-momentum from the HepMC3 particle
    const auto& mom = fp->momentum();

    // Convert each one from a HepMC3::GenParticle into a G4PrimaryParticle.
    // Do this by first setting the PDG code and the 4-momentum components.
    G4PrimaryParticle* particle;
    if (fp->pid() > 1000000000) {
    particle = new G4PrimaryParticle(
        fp->pid(),
        mom.px() * mom4_conv_factor,
        mom.py() * mom4_conv_factor,
        mom.pz() * mom4_conv_factor
    );
    }
    else{
      particle = new G4PrimaryParticle( fp->pid(),
      mom.px() * mom4_conv_factor,
      mom.py() * mom4_conv_factor,
      mom.pz() * mom4_conv_factor,
      mom.e() * mom4_conv_factor );
    }


    // Also set the charge of the G4PrimaryParticle appropriately.
    // Use the MARLEY utility function that checks for a stored charge
    // attribute before falling back to a PDG-code-based lookup.
    particle->SetCharge( marley_hepmc3::get_particle_charge( *fp ) );

    // Add the fully-initialized G4PrimaryParticle to the primary vertex
    vertex->SetPrimary( particle );
  }

  // The primary vertex has been fully populated with all final-state particles
  // from the MARLEY event. Add it to the G4Event object so that Geant4 can
  // begin tracking the particles through the simulated geometry.
  anEvent->AddPrimaryVertex( vertex );
}
G4ThreeVector MarleyPrimaryGeneratorAction::SamplePrimaryIronPosition()
{
  G4double L = 62.0 * m;
  G4double W = 15.1 * m;
  G4double H = 14.0 * m;
  G4double d = 1.2 * mm;

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

G4ThreeVector MarleyPrimaryGeneratorAction::SampleOuterIronPosition()
{

  G4double PrimaryIronThickness = 1.2 * mm;
  G4double PUThickness = 80 * cm;

  G4double L = 62.0 * m;
  G4double W = 15.1 * m;
  G4double H = 14.0 * m;
  G4double d = 54.63 * mm;

  //account for primary iron and PU Foam Layer
  L += 2*(PrimaryIronThickness + PUThickness);
  W += 2*(PrimaryIronThickness + PUThickness);
  H += 2*(PrimaryIronThickness + PUThickness);


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

std::array<double, 3> MarleyPrimaryGeneratorAction::SampleIsotropicDirection()
{
  double u = G4UniformRand();
  double v = G4UniformRand();

  double phi = 2*pi*u;
  double cos_theta = 2*v - 1;
  double sin_theta = std::sqrt(1 - cos_theta * cos_theta);

  double x = sin_theta * std::cos(phi);
  double y = sin_theta * std::sin(phi);
  double z = cos_theta;
  std::array<double, 3> nuDir = {x, y, z};
  return nuDir;
}

}