#ifndef EventAction_h
#define EventAction_h 1

#include "G4Types.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>
#include <map>

/// Event action class

class EventAction : public G4UserEventAction
{
public:
  EventAction() = default;
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event *event) override;
  void EndOfEventAction(const G4Event *event) override;
  G4int GetEventNo();

  // 添加边界穿越数据
  void AddBoundaryCrossing(G4int type, const G4ThreeVector &position);

private:
  // 存储当前事件的边界穿越数据：类型 -> 位置
  std::map<G4int, G4ThreeVector> fBoundaryCrossings;
  G4int eventCount = 0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
