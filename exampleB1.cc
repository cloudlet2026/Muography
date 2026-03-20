#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char **argv)
{
  // 如果是交互模式，则实例化 G4UIExecutive
  G4UIExecutive *ui = nullptr;
  if (argc == 1)
  {
    ui = new G4UIExecutive(argc, argv);
  }
  auto runManager = G4RunManagerFactory::CreateRunManager();
  // Force single-threaded mode to avoid allocator issues seen in MT runs
  // runManager->SetNumberOfThreads(1);
  G4int seed = 123;
  if (argc > 2)
    seed = atoi(argv[argc - 1]);

  // 选择随机引擎并设置种子
  // G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(seed);

  // 探测器构造
  auto detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);

  // 物理列表
  G4VModularPhysicsList *physicsList = new FTFP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  auto opticalPhysics = new G4OpticalPhysics();

  auto opticalParams = G4OpticalParameters::Instance();
  opticalParams->SetBoundaryInvokeSD(true);
  physicsList->RegisterPhysics(opticalPhysics);
  runManager->SetUserInitialization(physicsList);
  physicsList->DumpList();

  // 用户操作初始化
  runManager->SetUserInitialization(new ActionInitialization(detector));

  // 初始化可视化，仅用于交互式会话；在批处理模式下跳过以简化关闭
  G4VisManager *visManager = nullptr;
  if (ui)
  {
    visManager = new G4VisExecutive;
    visManager->Initialize();
  }

  // 获取用户界面管理器的指针
  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  // UImanager->ApplyCommand("/geant4e/limits/stepLength 0.1 mm");
  if (ui)
  {
    // 为交互模式定义（图形）用户界面终端
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
  }
  else
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  // 会话终止
  delete visManager;
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
