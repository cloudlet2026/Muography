# 缪子成像模拟项目 (Muography Simulation)

基于Geant4的缪子成像蒙特卡洛模拟与MLEM重建系统。

## 项目结构

```
Muography/
├── src/                    # Geant4源代码
│   ├── DetectorConstruction.cc   # 探测器几何定义
│   ├── PrimaryGeneratorAction.cc # 缪子源定义
│   ├── EventAction.cc            # 事件处理
│   ├── SteppingAction.cc         # 步进处理
│   └── RunAction.cc              # 运行处理
├── include/                # 头文件
├── muography_mlem.py       # MLEM重建脚本
├── muography_poca.py       # POCA成像脚本
├── analyze_output.C        # ROOT数据分析脚本
├── poca_imaging.C          # ROOT POCA成像脚本
├── CMakeLists.txt          # CMake配置
├── requirements.txt        # Python依赖
└── run1.mac, run2.mac      # Geant4宏文件
```

## 环境依赖

### 1. Geant4
- 版本：>= 11.0
- 需要UI和可视化支持

### 2. ROOT
- 版本：>= 6.20
- 用于数据输出

### 3. Python
- 版本：>= 3.9

```bash
pip install -r requirements.txt
```

## 编译步骤

```bash
cd build
cmake ..
make -j
```

## 运行模拟

### 方式1：使用宏文件

```bash
cd build
./exampleB1 run1.mac
```

### 方式2：交互模式

```bash
cd build
./exampleB1
```

### 输出文件

模拟完成后会在 `build/` 目录生成：
- `output.csv` - 探测器击中数据

## 数据处理与成像

### 1. MLEM重建 (Python)

```bash
python muography_mlem.py --input build/output.csv [选项]
```

#### 主要参数

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--input` | build/output.csv | 输入文件路径 |
| `--voxel` | 10.0 | 体素大小 (mm) |
| `--iterations` | 10 | EM迭代次数 |
| `--max-events` | 50000 | 最大事件数 |
| `--angle-min` | 1.0 | 最小散射角 (mrad) |
| `--angle-max` | 200.0 | 最大散射角 (mrad) |
| `--momentum` | 3000.0 | 缪子动量 (MeV/c) |
| `--radlen` | 304000.0 | 参考辐射长度 (mm) |
| `--no-plot` | - | 跳过绑图 |

#### 示例

```bash
# 基本运行
python muography_mlem.py

# 自定义参数
python muography_mlem.py --voxel 5.0 --iterations 20 --max-events 100000
```

### 2. POCA成像 (Python)

```bash
python muography_poca.py --input build/output.csv
```

### 3. ROOT数据分析 (analyze_output.C)

使用ROOT进行POCA分析和统计：

```bash
# 方式1：直接执行脚本
root -l -q 'analyze_output.C("build/output.csv")'

# 方式2：使用默认路径
root -l -q analyze_output.C
```

#### 输出内容

ROOT分析会输出以下统计信息：
- 总事件数和有效事件数
- POCA距离统计（平均值、最小值、最大值）
- 散射角统计（mrad）
- POCA点空间分布范围

## 输出结果

输出会在 `picture/` 目录生成：

## 几何配置

当前探测器配置（在 `src/DetectorConstruction.cc` 中定义）：

```
Z轴布局 (mm):
  800 ────────  模块1 (L1)
  600 ────────  模块2 (L2)
       │
       │  ROI区域
       │    ┌────────┐
       │    │  铅罐  │  (高150mm, 外径100mm, 壁厚20mm)
       │    └────────┘	
       │
 -600 ────────  模块3 (L3)
 -800 ────────  模块4 (L4)
```

