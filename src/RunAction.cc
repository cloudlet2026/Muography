#include "RunAction.hh"
#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "SteppingAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction *eventaction)
    : fEventAction(eventaction)
{
  // auto analysisManager = G4AnalysisManager::Instance();
  // analysisManager->SetDefaultFileType("root");
  // // If the filename extension is not provided, the default file type (root)
  // // will be used for all files specified without extension.
  // analysisManager->SetVerboseLevel(1);
  // analysisManager->SetNtupleMerging(true);
  // // Note: merging ntuples is available only with Root output
  // analysisManager->SetFileName("pocahisto");
  // // Book histograms, ntuple
  // if (fEventAction)
  // {
  //   // Creating ntuple
  //   analysisManager->CreateNtuple("ntuple", "index of detector"); // 创建元组
  // }
  // Set ntuple output file
  // analysisManager->SetNtupleFileName(0, "psfntuple");
  // 若不存在则创建
  if (!std::ifstream("output.txt"))
  {
    std::ofstream outFile("output.txt");
    outFile.close();
  }
  else // 若存在则清空
  {
    std::ofstream outFile("output.txt", std::ios::trunc);
    // 清空文件内容
    outFile << "";
    outFile.close();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run *)
{

  // auto analysisManager = G4AnalysisManager::Instance();
  // analysisManager->Reset();
  // analysisManager->OpenFile(); // do NOT insert fileName, otherwise it will not accept the name from the .mac file
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run *)
{

  // auto analysisManager = G4AnalysisManager::Instance();
  // analysisManager->Write();
  // analysisManager->CloseFile(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
