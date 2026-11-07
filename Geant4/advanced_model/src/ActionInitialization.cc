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
/// \file ActionInitialization.cc
/// \brief Implementation of the B1::ActionInitialization class

#include "ActionInitialization.hh"
#include <cstdlib>
#include <string>
#include "G4Exception.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "MarleyPrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  auto runAction = new RunAction;
  SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Set the environment variable FDSTEEL_CONFIG in bash
// = IBD for only Geant4
void ActionInitialization::Build() const
{

  const char* config = std::getenv("FDSTEEL_CONFIG");

  if (!config)
  {
      G4Exception(
          "ActionInitialization::Build",
          "MissingConfig",
          FatalException,
          "FDSTEEL_CONFIG is not set!"
      );
      return;
  }

  if (std::string(config) == "IBD")
  {
      SetUserAction(new PrimaryGeneratorAction);
  }
  else
  {
      SetUserAction(new MarleyPrimaryGeneratorAction(config));
  }

  SetUserAction(new RunAction);

  auto eventAction = new EventAction;
  SetUserAction(eventAction);

  SetUserAction(new SteppingAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
