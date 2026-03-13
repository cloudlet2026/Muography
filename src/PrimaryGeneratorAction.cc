#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4UImanager.hh"
#include "G4AutoLock.hh"

G4Mutex gen_mutex = G4MUTEX_INITIALIZER;
G4bool PrimaryGeneratorAction::fSeaLevel = true;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fParticleGun = new G4GeneralParticleSource();
  fParticleGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("mu-"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;

  if (fFlux)
  {
    delete fFlux;
    fFlux = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double getFlux(double muonEnergy, double theta)
{
  const double P1 = 0.102573;
  const double P2 = -0.068287;
  const double P3 = 0.958633;
  const double P4 = 0.0407253;
  const double P5 = 0.817285;

  double cosTheta = cos(theta);
  double cosThetaStar2 = (cosTheta * cosTheta + P1 * P1 + P2 * pow(cosTheta, P3) + P4 * pow(cosTheta, P5)) / (1. + P1 * P1 + P2 + P4);
  double cosThetaStar = sqrt(cosThetaStar2);

  double Emu = muonEnergy; // GeV

  double term1 = 1. / (1 + 1.1 * Emu * cosThetaStar / 115);
  double term2 = 0.054 / (1 + 1.1 * Emu * cosThetaStar / 850);

  double flux = 0.14 * pow(Emu * (1. + 3.64 / (Emu * pow(cosThetaStar, 1.29))), -2.7) * (term1 + term2);

  return flux * 10000; // 1/cm^2 -> 1/m^2
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

double funFlux(double *x, double *par)
{
  double muonEnergy = x[0];
  double theta = x[1] / 180. * TMath::Pi();

  return getFlux(muonEnergy, theta);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  if (!fFlux)
  {
    G4AutoLock l(&gen_mutex);
    if (!fFlux)
      fFlux = new TF2("f0", funFlux, 0.1, 2000, 0, 89);
  }
  double energy, theta, theta_rad;
  G4double phi = G4UniformRand() * 2.0 * CLHEP::pi;
  // 从二维分布中采样能量和天顶角theta
  if (fFlux)
    fFlux->GetRandom2(energy, theta);
  if (!fSeaLevel)
  {
    energy = 3.; // 固定能量
    theta = 0.;  // 垂直入射
  }
  // 天顶角theta定义为与z-轴的夹角，定义xOy平面为水平面
  theta_rad = theta / 180. * TMath::Pi();
  G4double z_dir = -cos(theta_rad);
  G4double y_dir = sin(theta_rad) * sin(phi);
  G4double x_dir = sin(theta_rad) * cos(phi);
  // 归一化方向矢量
  G4ThreeVector fParticleGun_direction(x_dir, y_dir, z_dir);
  fParticleGun_direction = fParticleGun_direction.unit();
  // 在探测器上方某个平面随机生成位置
  // 扩大源面覆盖范围，覆盖±400 mm 的探测区域
  G4double fParticleGun_size_x = 100. * cm;
  G4double fParticleGun_size_y = 100. * cm;
  G4double fParticleGun_size_z = 81. * cm;
  double x_pos = (G4UniformRand() - 0.5) * fParticleGun_size_x;
  double y_pos = (G4UniformRand() - 0.5) * fParticleGun_size_y;
  double z_pos = fParticleGun_size_z;
  G4ThreeVector fParticleGun_position(x_pos, y_pos, z_pos);
  auto current = fParticleGun->GetCurrentSource();

  if (!fFlux)
  {
    G4AutoLock l(&gen_mutex);
    if (!fFlux)
      fFlux = new TF2("f0", funFlux, 0.1, 2000, 0, 89);
  }
  if (current)
  {
    if (current->GetEneDist())
    {
      current->GetEneDist()->SetMonoEnergy(energy * GeV);
      current->GetEneDist()->SetEnergyDisType("Mono");
    }
    if (current->GetAngDist())
    {
      current->GetAngDist()->SetParticleMomentumDirection(G4ThreeVector(x_dir, y_dir, z_dir));
    }
    if (current->GetPosDist())
    {
      current->GetPosDist()->SetCentreCoords(G4ThreeVector(x_pos, y_pos, z_pos));
    }
  }

  // The code behind this line is not thread safe because polarization and time are randomly selected and GPS properties are global
  // 代码不是线程安全的，因为极化和时间是随机选择的，GPS属性是全局的

  G4AutoLock l(&gen_mutex);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
