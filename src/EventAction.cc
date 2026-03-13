#include "EventAction.hh"
#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *event)
{
  // 清空上一个事件的数据
  fBoundaryCrossings.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *evt)
{

  // 检查是否包含所有四种边界类型（1、2、3、4）
  if (fBoundaryCrossings.size() == 4 &&
      fBoundaryCrossings.count(1) &&
      fBoundaryCrossings.count(2) &&
      fBoundaryCrossings.count(3) &&
      fBoundaryCrossings.count(4))
  {
    eventCount++;
    // 按顺序输出四个边界
    for (G4int i = 1; i <= 4; i++)
    {
      const G4ThreeVector &pos = fBoundaryCrossings[i];
      // 写入output.txt 覆盖模式
      std::ofstream outFile("output.txt", std::ios::app);
      outFile << eventCount << "\t" << i << "\t" << pos << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddBoundaryCrossing(G4int type, const G4ThreeVector &position)
{
  // 存储边界穿越数据（如果同一类型多次穿越，只保留第一次）
  if (fBoundaryCrossings.count(type) == 0)
  {
    fBoundaryCrossings[type] = position;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int EventAction::GetEventNo()
{
  return G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
}
