#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

namespace B1
{

RunAction::RunAction()
{
  // Print progress every 1000 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);

  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->SetVerboseLevel(1);

  // Needed for ROOT output in multi-threaded mode
  analysisManager->SetNtupleMerging(true);

  // Output file
  analysisManager->SetFileName("FDSteel.root");

  // -------------------------
  // Histogram
  // -------------------------

  analysisManager->CreateH1(
      "EdepLAr",
      "Energy deposited in liquid argon",
      2000,
      0.0,
      200.0 * MeV
  );

  // -------------------------
  // Ntuple
  // -------------------------

  analysisManager->CreateNtuple(
      "FDSteel",
      "FDSteel event data"
  );

  analysisManager->CreateNtupleDColumn("EdepLAr");

  analysisManager->FinishNtuple();
}


void RunAction::BeginOfRunAction(const G4Run*)
{
  auto analysisManager = G4AnalysisManager::Instance();

  analysisManager->OpenFile();
}


void RunAction::EndOfRunAction(const G4Run*)
{
  auto analysisManager = G4AnalysisManager::Instance();

  if (analysisManager->GetH1(0))
  {
    if (IsMaster())
    {
      G4cout << G4endl;
      G4cout << "========== FDSteel Run Summary ==========" << G4endl;

      G4cout
          << "EdepLAr mean = "
          << G4BestUnit(
              analysisManager->GetH1(0)->mean(),
              "Energy"
          )
          << G4endl;

      G4cout
          << "EdepLAr RMS  = "
          << G4BestUnit(
              analysisManager->GetH1(0)->rms(),
              "Energy"
          )
          << G4endl;

      G4cout << "=========================================" << G4endl;
    }
  }

  analysisManager->Write();
  analysisManager->CloseFile();
}

}  // namespace B1