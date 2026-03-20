#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph2D.h>
#include <TStyle.h>

struct EventData {
    Int_t eventId;
    Int_t moduleId;
    Double_t x;
    Double_t y;
};

struct Point {
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

void analyze_output(const char* inputFile = "build/output.csv") {
    // 打开输入文件
    std::ifstream inFile(inputFile);
    if (!inFile.is_open()) {
        std::cerr << "Error: Cannot open file " << inputFile << std::endl;
        return;
    }
    
    // 读取表头
    std::string header;
    std::getline(inFile, header);
    
    // 存储事件数据
    std::map<Int_t, std::map<Int_t, EventData>> eventData;
    
    // 解析数据
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
    
    // 几何常数（毫米）- 与仿真一致
    const Double_t z1 = 800.0;
    const Double_t z2 = 600.0;
    const Double_t z3 = -600.0;
    const Double_t z4 = -800.0;
    
    // 存储POCA点和相关指标
    Int_t totalEvents = 0;
    Int_t validEvents = 0;
    std::vector<Point> pocaPoints;
    std::vector<Double_t> pocaDistances;
    std::vector<Double_t> scatteringAngles;
    
    // 处理每个事件
    for (auto& eventEntry : eventData) {
        Int_t eventId = eventEntry.first;
        auto& modules = eventEntry.second;
        totalEvents++;
        
        // 检查是否有所有4个层的数据
        if (modules.find(1) == modules.end() || modules.find(2) == modules.end() ||
            modules.find(3) == modules.end() || modules.find(4) == modules.end()) {
            continue;
        }
        
        // 获取各层数据
        EventData& m1 = modules[1];
        EventData& m2 = modules[2];
        EventData& m3 = modules[3];
        EventData& m4 = modules[4];
        
        // 入射直线：模块1 -> 模块2 (从上往下)
        Line incidentLine = {
            {m1.x, m1.y, z1},
            {m2.x, m2.y, z2}
        };
        
        // 出射直线：模块3 -> 模块4 (从上往下)
        Line exitLine = {
            {m3.x, m3.y, z3},
            {m4.x, m4.y, z4}
        };
        
        // 计算POCA距离
        Double_t distance = incidentLine.distanceToLine(exitLine);
        
        // 计算散射角
        Point dir1 = incidentLine.direction();
        Point dir2 = exitLine.direction();
        
        Double_t dot = dir1.x * dir2.x + dir1.y * dir2.y + dir1.z * dir2.z;
        Double_t mag1 = std::sqrt(dir1.x * dir1.x + dir1.y * dir1.y + dir1.z * dir1.z);
        Double_t mag2 = std::sqrt(dir2.x * dir2.x + dir2.y * dir2.y + dir2.z * dir2.z);
        
        if (mag1 < 1e-10 || mag2 < 1e-10) continue;
        
        Double_t cosAngle = dot / (mag1 * mag2);
        if (cosAngle > 1.0) cosAngle = 1.0;
        if (cosAngle < -1.0) cosAngle = -1.0;
        Double_t scatteringAngle = std::acos(cosAngle);
        
        // 计算POCA点
        Point poca = incidentLine.closestPointToLine(exitLine);
        
        // 存储数据
        validEvents++;
        pocaPoints.push_back(poca);
        pocaDistances.push_back(distance);
        scatteringAngles.push_back(scatteringAngle);
    }
    
    // 计算统计指标
    Double_t totalPocaDistance = 0;
    Double_t totalScatteringAngle = 0;
    Double_t minPocaDistance = 1e9;
    Double_t maxPocaDistance = 0;
    Double_t minScatteringAngle = 1e9;
    Double_t maxScatteringAngle = 0;
    
    for (size_t i = 0; i < pocaDistances.size(); i++) {
        totalPocaDistance += pocaDistances[i];
        totalScatteringAngle += scatteringAngles[i];
        
        if (pocaDistances[i] < minPocaDistance) minPocaDistance = pocaDistances[i];
        if (pocaDistances[i] > maxPocaDistance) maxPocaDistance = pocaDistances[i];
        if (scatteringAngles[i] < minScatteringAngle) minScatteringAngle = scatteringAngles[i];
        if (scatteringAngles[i] > maxScatteringAngle) maxScatteringAngle = scatteringAngles[i];
    }
    
    Double_t avgPocaDistance = (validEvents > 0) ? totalPocaDistance / validEvents : 0;
    Double_t avgScatteringAngle = (validEvents > 0) ? totalScatteringAngle / validEvents : 0;
    
    // 计算POCA点分布范围
    Double_t minX = 1e9, maxX = -1e9;
    Double_t minY = 1e9, maxY = -1e9;
    Double_t minZ = 1e9, maxZ = -1e9;
    
    for (const auto& p : pocaPoints) {
        if (p.x < minX) minX = p.x;
        if (p.x > maxX) maxX = p.x;
        if (p.y < minY) minY = p.y;
        if (p.y > maxY) maxY = p.y;
        if (p.z < minZ) minZ = p.z;
        if (p.z > maxZ) maxZ = p.z;
    }
    
    // 打印指标信息
    std::cout << "=== 成像数据分析结果 ===" << std::endl;
    std::cout << "输入文件: " << inputFile << std::endl;
    std::cout << "总事件数: " << totalEvents << std::endl;
    std::cout << "有效事件数: " << validEvents << std::endl;
    std::cout << "有效事件比例: " << (totalEvents > 0 ? (Double_t)validEvents / totalEvents * 100 : 0) << "%" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "=== POCA 距离统计 ===" << std::endl;
    std::cout << "平均距离: " << avgPocaDistance << " mm" << std::endl;
    std::cout << "最小距离: " << minPocaDistance << " mm" << std::endl;
    std::cout << "最大距离: " << maxPocaDistance << " mm" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "=== 散射角统计 ===" << std::endl;
    std::cout << "平均散射角: " << avgScatteringAngle * 1000 << " mrad" << std::endl;
    std::cout << "最小散射角: " << minScatteringAngle * 1000 << " mrad" << std::endl;
    std::cout << "最大散射角: " << maxScatteringAngle * 1000 << " mrad" << std::endl;
    std::cout << "" << std::endl;
    std::cout << "=== POCA 点分布范围 ===" << std::endl;
    std::cout << "X范围: [" << minX << ", " << maxX << "] mm" << std::endl;
    std::cout << "Y范围: [" << minY << ", " << maxY << "] mm" << std::endl;
    std::cout << "Z范围: [" << minZ << ", " << maxZ << "] mm" << std::endl;
    std::cout << "" << std::endl;
    
}
