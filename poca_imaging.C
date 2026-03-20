#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <TFile.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>

struct Point {
    Double_t x, y, z;
};

struct EventData {
    Int_t eventId;
    Int_t moduleId;
    Double_t x, y, z;
};

struct Line {
    Point p1, p2;
    
    Point direction() const {
        return {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
    }
    
    Point closestPointToLine(const Line& other) const {
        Point d1 = direction();
        Point d2 = other.direction();
        
        Double_t a = d1.x * d1.x + d1.y * d1.y + d1.z * d1.z;
        Double_t b = d1.x * d2.x + d1.y * d2.y + d1.z * d2.z;
        Double_t c = d2.x * d2.x + d2.y * d2.y + d2.z * d2.z;
        Double_t d = d1.x * (p1.x - other.p1.x) + d1.y * (p1.y - other.p1.y) + d1.z * (p1.z - other.p1.z);
        Double_t e = d2.x * (p1.x - other.p1.x) + d2.y * (p1.y - other.p1.y) + d2.z * (p1.z - other.p1.z);
        
        Double_t denom = a * c - b * b;
        
        if (std::abs(denom) < 1e-10) return p1;
        
        Double_t sc = (b * d - a * e) / denom;
        Double_t tc = (c * d - b * e) / denom;
        
        Point closest1 = {p1.x + sc * d1.x, p1.y + sc * d1.y, p1.z + sc * d1.z};
        Point closest2 = {other.p1.x + tc * d2.x, other.p1.y + tc * d2.y, other.p1.z + tc * d2.z};
        
        return {(closest1.x + closest2.x) / 2, (closest1.y + closest2.y) / 2, (closest1.z + closest2.z) / 2};
    }
    
    Double_t distanceToLine(const Line& other) const {
        Point d1 = direction();
        Point d2 = other.direction();
        
        Point r = {p1.x - other.p1.x, p1.y - other.p1.y, p1.z - other.p1.z};
        
        Double_t a = d1.x * d1.x + d1.y * d1.y + d1.z * d1.z;
        Double_t b = d1.x * d2.x + d1.y * d2.y + d1.z * d2.z;
        Double_t c = d2.x * d2.x + d2.y * d2.y + d2.z * d2.z;
        Double_t d = d1.x * r.x + d1.y * r.y + d1.z * r.z;
        Double_t e = d2.x * r.x + d2.y * r.y + d2.z * r.z;
        
        Double_t denom = a * c - b * b;
        
        if (std::abs(denom) < 1e-10) {
            return std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
        }
        
        Double_t sc = (b * d - a * e) / denom;
        Double_t tc = (c * d - b * e) / denom;
        
        Point closest1 = {p1.x + sc * d1.x, p1.y + sc * d1.y, p1.z + sc * d1.z};
        Point closest2 = {other.p1.x + tc * d2.x, other.p1.y + tc * d2.y, other.p1.z + tc * d2.z};
        
        Double_t dx = closest1.x - closest2.x;
        Double_t dy = closest1.y - closest2.y;
        Double_t dz = closest1.z - closest2.z;
        
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }
};

Double_t calculateScatteringAngle(const Line& incidentLine, const Line& exitLine) {
    Point dir1 = incidentLine.direction();
    Point dir2 = exitLine.direction();
    
    Double_t dot = dir1.x * dir2.x + dir1.y * dir2.y + dir1.z * dir2.z;
    Double_t mag1 = std::sqrt(dir1.x * dir1.x + dir1.y * dir1.y + dir1.z * dir1.z);
    Double_t mag2 = std::sqrt(dir2.x * dir2.x + dir2.y * dir2.y + dir2.z * dir2.z);
    
    if (mag1 < 1e-10 || mag2 < 1e-10) return 0;
    
    Double_t cosAngle = dot / (mag1 * mag2);
    if (cosAngle > 1.0) cosAngle = 1.0;
    if (cosAngle < -1.0) cosAngle = -1.0;
    
    return std::acos(cosAngle);
}

Double_t calculateWeight(Double_t distance, Double_t scatteringAngle) {
    return scatteringAngle * 1000.0;
}

void poca_imaging(const char* inputFile = "build/output.csv",
                      const char* outputFile = "poca_image.root") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetPalette(kRainBow);
    
    std::cout << "=== POCA 图像重建 ===" << std::endl;
    std::cout << "输入文件: " << inputFile << std::endl;
    
    std::ifstream inFile(inputFile);
    if (!inFile.is_open()) {
        std::cerr << "Error: Cannot open file " << inputFile << std::endl;
        return;
    }
    
    std::string header;
    std::getline(inFile, header);
    
    std::map<Int_t, std::map<Int_t, EventData>> eventData;
    
    std::string line;
    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        EventData data;
        char comma;
        ss >> data.eventId >> comma >> data.moduleId >> comma 
           >> data.x >> comma >> data.y;
        eventData[data.eventId][data.moduleId] = data;
    }
    inFile.close();
    
    const Double_t z1 = 800.0;
    const Double_t z2 = 600.0;
    const Double_t z3 = -600.0;
    const Double_t z4 = -800.0;
    
    const Double_t minX = -300.0, maxX = 300.0;
    const Double_t minY = -300.0, maxY = 300.0;
    const Double_t minZ = -200.0, maxZ = 200.0;
    const Int_t bins = 60;
    
    Int_t totalEvents = 0;
    Int_t validEvents = 0;
    std::vector<Point> pocaPoints;
    std::vector<Double_t> pocaDistances;
    std::vector<Double_t> scatteringAngles;
    std::vector<Double_t> weights;
    
    for (auto& eventEntry : eventData) {
        totalEvents++;
        auto& modules = eventEntry.second;
        
        if (modules.find(1) == modules.end() || modules.find(2) == modules.end() ||
            modules.find(3) == modules.end() || modules.find(4) == modules.end()) {
            continue;
        }
        
        EventData& m1 = modules[1];
        EventData& m2 = modules[2];
        EventData& m3 = modules[3];
        EventData& m4 = modules[4];
        
        Line incidentLine = {{m1.x, m1.y, z1}, {m2.x, m2.y, z2}};
        Line exitLine = {{m4.x, m4.y, z4}, {m3.x, m3.y, z3}};
        Line exitLineReal = {{m3.x, m3.y, z3}, {m4.x, m4.y, z4}};
        
        Point poca = incidentLine.closestPointToLine(exitLine);
        Double_t distance = incidentLine.distanceToLine(exitLine);
        Double_t scatteringAngle = calculateScatteringAngle(incidentLine, exitLineReal);
        Double_t weight = calculateWeight(distance, scatteringAngle);
        
        if (weight > 0) {
            validEvents++;
            pocaPoints.push_back(poca);
            pocaDistances.push_back(distance);
            scatteringAngles.push_back(scatteringAngle * 1000.0);
            weights.push_back(weight);
        }
    }
    
    std::cout << "总事件: " << totalEvents << std::endl;
    std::cout << "有效事件: " << validEvents << std::endl;
    
    TFile* outFile = new TFile(outputFile, "RECREATE");
    
    TH3D* hPOCA3D = new TH3D("hPOCA3D", "POCA 3D;X (mm);Y (mm);Z (mm)",
                              bins, minX, maxX, bins, minY, maxY, bins, minZ, maxZ);
    TH2D* hXY = new TH2D("hXY", "POCA XY;X (mm);Y (mm)", bins, minX, maxX, bins, minY, maxY);
    TH2D* hXZ = new TH2D("hXZ", "POCA XZ;X (mm);Z (mm)", bins, minX, maxX, bins, minZ, maxZ);
    TH2D* hYZ = new TH2D("hYZ", "POCA YZ;Y (mm);Z (mm)", bins, minY, maxY, bins, minZ, maxZ);
    
    TH1D* hDist = new TH1D("hDist", "POCA Distance;Distance (mm);Counts", 100, 0, 10);
    TH1D* hAngle = new TH1D("hAngle", "Scattering Angle;Angle (mrad);Counts", 100, 0, 50);
    
    for (size_t i = 0; i < pocaPoints.size(); i++) {
        Point& p = pocaPoints[i];
        if (p.x >= minX && p.x <= maxX && p.y >= minY && p.y <= maxY && 
            p.z >= minZ && p.z <= maxZ) {
            hPOCA3D->Fill(p.x, p.y, p.z, weights[i]);
            hXY->Fill(p.x, p.y, weights[i]);
            hXZ->Fill(p.x, p.z, weights[i]);
            hYZ->Fill(p.y, p.z, weights[i]);
        }
        hDist->Fill(pocaDistances[i]);
        hAngle->Fill(scatteringAngles[i]);
    }
    
    TCanvas* c3D = new TCanvas("c3D", "POCA 3D", 800, 600);
    hPOCA3D->Draw("BOX");
    c3D->Write();
    
    TCanvas* cXY = new TCanvas("cXY", "POCA XY", 800, 600);
    hXY->Draw("COLZ");
    cXY->Write();
    
    TCanvas* cXZ = new TCanvas("cXZ", "POCA XZ", 800, 600);
    hXZ->Draw("COLZ");
    cXZ->Write();
    
    TCanvas* cYZ = new TCanvas("cYZ", "POCA YZ", 800, 600);
    hYZ->Draw("COLZ");
    cYZ->Write();
    
    TCanvas* cDist = new TCanvas("cDist", "POCA Distance", 800, 600);
    hDist->Draw();
    cDist->Write();
    
    TCanvas* cAngle = new TCanvas("cAngle", "Scattering Angle", 800, 600);
    hAngle->Draw();
    cAngle->Write();
    
    outFile->Write();
    outFile->Close();
    
    std::cout << "输出文件: " << outputFile << std::endl;
}
