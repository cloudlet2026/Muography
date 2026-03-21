#!/usr/bin/env python3

from __future__ import annotations
import argparse
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import font_manager
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

# ----------------几何常数 ----------------
PLANE_Z = {1: 800.0, 2: 600.0, 3: -600.0, 4: -800.0}  # mm
# Z 多深度切片（指定 Z=-40,-20,-5,5,20,40 cm），2x3 布局
target_z_vals = [-100, -30, -10, 0, 10, 30]
# ROI 范围（可以通过 CLI 覆盖）
DEFAULT_X_RANGE = (-600,600)
DEFAULT_Y_RANGE = (-600,600)
DEFAULT_Z_RANGE = (PLANE_Z[3], PLANE_Z[2])

# ----------------高地参数 ----------------
HIGHLAND_K = 13.6  # MeV
DEFAULT_MUON_MOMENTUM = 3000.0  # MeV/c，典型的宇宙μ子
# 参考辐射长度：使用空气作为基准（约304m = 304000mm）
# 密度值表示相对于空气的散射能力倍数
# 例如：密度=100表示该体素散射能力是空气的100倍
DEFAULT_RADIATION_LENGTH = 304000.0  # mm, air
font_manager.fontManager.addfont("fonts/TimesSimSunRegular.ttf")
plt.rcParams["font.family"] = "TimesSimSun"  # 支持中字体

# ---------------- 几何参数 ----------------
Z1 = 800.0   # 模块1
Z2 = 600.0    # 模块2
Z3 = -600.0   # 模块3
Z4 = -800.0  # 模块4
# ------------------------------------------

def parse_output(path: Path, max_events: int | None = None) -> Dict[int, Dict[int, np.ndarray]]:
    """将output.csv解析为{event_id: {layer_id: np.array([x,y,z])}}。
    如果提供，则在 max_events（完成事件）之后停止。
    """
    events: Dict[int, Dict[int, np.ndarray]] = {}
    current_id = None
    for line in path.open():
        line = line.strip()
        if not line:
            continue
        # 跳过CSV表头
        if line.startswith('eventId'):
            continue
        parts = line.split(',')
        if len(parts) < 4:
            continue
        event_id = int(parts[0])
        layer_id = int(parts[1])
        x = float(parts[2])
        y = float(parts[3])
        
        # 根据layer_id设置z坐标
        if layer_id == 1:
            z = Z1
        elif layer_id == 2:
            z = Z2
        elif layer_id == 3:
            z = Z3
        elif layer_id == 4:
            z = Z4
        else:
            z = 0.0

        if event_id not in events:
            events[event_id] = {}
        events[event_id][layer_id] = np.array([x, y, z], dtype=float)

        # 一旦收集到足够的完整事件，就提前停止
        if max_events and len(events) >= max_events:
            # 仅当当前事件完成时才停止
            if len(events[event_id]) == 4:
                break
    return events

# ----------------数学求解器 ----------------

def scattering_angle_mrad(dir_top: np.ndarray, dir_bottom_down: np.ndarray) -> float:
    """两个方向之间的夹角，以 mrad 为单位。"""
    u1 = dir_top / (np.linalg.norm(dir_top) + 1e-12)
    u2 = dir_bottom_down / (np.linalg.norm(dir_bottom_down) + 1e-12)
    cos_theta = float(np.clip(np.dot(u1, u2), -1.0, 1.0))
    theta = math.acos(cos_theta)
    return theta * 1000.0

def highland_pdf(angle_mrad: float, path_length_mm: float,
                 momentum_mev: float, radiation_length_mm: float) -> float:
    """获取 高地公式 的高斯概率密度（角度以 mrad 为单位）。"""
    x_over_x0 = max(path_length_mm / max(radiation_length_mm, 1e-6), 1e-9)
    theta0_rad = (HIGHLAND_K / momentum_mev) * math.sqrt(x_over_x0)
    theta0_rad *= 1.0 + 0.038 * math.log(x_over_x0)
    theta0_mrad = max(theta0_rad * 1000.0, 1e-9)
    exponent = -0.5 * (angle_mrad / theta0_mrad) ** 2
    norm = theta0_mrad * math.sqrt(2.0 * math.pi)
    return math.exp(exponent) / norm

# ---------------- 网格和追踪 ----------------

class VoxelGrid:
    def __init__(self, x_range: Tuple[float, float], y_range: Tuple[float, float],
                 z_range: Tuple[float, float], voxel_size: float):
        self.x_range = x_range
        self.y_range = y_range
        self.z_range = z_range
        self.voxel_size = voxel_size
        self.origin = np.array([x_range[0], y_range[0], z_range[0]], dtype=float)
        self.nx = int(math.ceil((x_range[1] - x_range[0]) / voxel_size))
        self.ny = int(math.ceil((y_range[1] - y_range[0]) / voxel_size))
        self.nz = int(math.ceil((z_range[1] - z_range[0]) / voxel_size))
        self.density = np.ones((self.nx, self.ny, self.nz), dtype=float)

    @property
    def bounds_min(self) -> np.ndarray:
        return self.origin

    @property
    def bounds_max(self) -> np.ndarray:
        return self.origin + np.array([self.nx, self.ny, self.nz], dtype=float) * self.voxel_size

    def inside(self, point: np.ndarray) -> bool:
        x, y, z = point
        return (self.x_range[0] <= x <= self.x_range[1] and
                self.y_range[0] <= y <= self.y_range[1] and
                self.z_range[0] <= z <= self.z_range[1])

    def index(self, point: np.ndarray) -> Tuple[int, int, int]:
        rel = (point - self.origin) / self.voxel_size
        return int(rel[0]), int(rel[1]), int(rel[2])

    def flat_index(self, idx3d: Tuple[int, int, int]) -> int:
        ix, iy, iz = idx3d
        return (ix * self.ny + iy) * self.nz + iz

    def total_voxels(self) -> int:
        return self.nx * self.ny * self.nz


def _clip_line_to_box(p: np.ndarray, d: np.ndarray, bmin: np.ndarray, bmax: np.ndarray) -> Tuple[float, float]:
    # """Return (t_enter, t_exit) for line p + t d clipped to AABB."""
    t_enter, t_exit = -math.inf, math.inf
    for i in range(3):
        if abs(d[i]) < 1e-12:
            if p[i] < bmin[i] or p[i] > bmax[i]:
                return math.inf, -math.inf
            continue
        t1 = (bmin[i] - p[i]) / d[i]
        t2 = (bmax[i] - p[i]) / d[i]
        t_min = min(t1, t2)
        t_max = max(t1, t2)
        t_enter = max(t_enter, t_min)
        t_exit = min(t_exit, t_max)
        if t_exit < t_enter:
            return math.inf, -math.inf
    return t_enter, t_exit


def trace_line(p: np.ndarray, d: np.ndarray, grid: VoxelGrid) -> Dict[int, float]:
    # """Siddon-style traversal of line p + t d through the ROI box.
    # Returns {flat_index: path_length_mm} for the segment inside ROI.
    # """
    if np.linalg.norm(d) < 1e-9:
        return {}

    bmin, bmax = grid.bounds_min, grid.bounds_max
    t_enter, t_exit = _clip_line_to_box(p, d, bmin, bmax)
    if not math.isfinite(t_enter) or t_exit <= t_enter:
        return {}

    # 限制为 ROI 内的段
    start = p + t_enter * d
    end = p + t_exit * d
    seg_dir = end - start
    seg_len = float(np.linalg.norm(seg_dir))
    if seg_len < 1e-6:
        return {}
    seg_dir_unit = seg_dir / seg_len

    # 初始体素索引
    idx = np.floor((start - grid.origin) / grid.voxel_size).astype(int)
    idx = np.clip(idx, 0, [grid.nx - 1, grid.ny - 1, grid.nz - 1])

    step = np.sign(seg_dir_unit).astype(int)
    step[step == 0] = 1

    # 每个轴的下一个边界（世界坐标）
    next_boundary = grid.origin + (idx + (step > 0)) * grid.voxel_size
    t_max = (next_boundary - start) / (seg_dir_unit + 1e-12)  # mm along line
    t_delta = grid.voxel_size / (np.abs(seg_dir_unit) + 1e-12)

    lengths: Dict[int, float] = {}
    t = 0.0
    max_voxels = grid.total_voxels() * 2  # safety

    for _ in range(max_voxels):
        if t >= seg_len:
            break
        axis = int(np.argmin(t_max))
        t_next = float(t_max[axis])
        path = min(t_next, seg_len) - t
        flat = grid.flat_index((int(idx[0]), int(idx[1]), int(idx[2])))
        lengths[flat] = lengths.get(flat, 0.0) + path
        t = t_next
        idx[axis] += step[axis]
        if (idx[0] < 0 or idx[0] >= grid.nx or
                idx[1] < 0 or idx[1] >= grid.ny or
                idx[2] < 0 or idx[2] >= grid.nz):
            break
        t_max[axis] += t_delta[axis]

    return lengths

# ---------------- EM-ML ----------------

class EventData:
    def __init__(self, contributions: Dict[int, float], angle_mrad: float, total_length: float):
        self.contributions = contributions
        self.angle_mrad = angle_mrad
        self.total_length = total_length


class EMLReconstructor:
    def __init__(self, grid: VoxelGrid, events: List[EventData],
                 momentum_mev: float, radiation_length_mm: float):
        self.grid = grid
        self.events = events
        self.momentum = momentum_mev
        self.rad_len = radiation_length_mm

    def compute_effective_path(self, evt: EventData, density_flat: np.ndarray) -> float:
        """
        计算体素密度加权的等效路径长度。
        等效路径 = Σ(路径长度_i × 密度_i)
        这个值用于计算等效的 x/X0
        """
        effective_path = 0.0
        for vidx, path_len in evt.contributions.items():
            effective_path += path_len * density_flat[vidx]
        return effective_path

    def iterate(self) -> None:
        n_vox = self.grid.total_voxels()
        density_flat = self.grid.density.reshape(-1)
        new_density = np.zeros_like(density_flat)

        for evt in self.events:
            if not evt.contributions:
                continue
            
            effective_path = self.compute_effective_path(evt, density_flat)
            
            if effective_path <= 1e-12:
                continue
            
            weight_angle = highland_pdf(evt.angle_mrad, effective_path,
                                        self.momentum, self.rad_len)
            
            expected = effective_path
            norm = weight_angle / expected
            for vidx, path_len in evt.contributions.items():
                new_density[vidx] += density_flat[vidx] * path_len * norm

        mean_val = float(new_density.mean())
        if mean_val > 0:
            new_density /= mean_val
        self.grid.density = new_density.reshape(self.grid.density.shape)

    def run(self, iterations: int) -> None:
        for i in range(iterations):
            self.iterate()
            print(f"迭代 {i + 1}/{iterations} 完成。平均值={self.grid.density.mean():.3f} 最大值={self.grid.density.max():.3f}")

# ----------------流水线 ----------------

def save_density_plots(grid: VoxelGrid, prefix: str = "./picture/mlem",
                       density_percentile: float = 90.0,
                       dbscan_eps_factor: float = 1.0,
                       dbscan_min_samples: int = 3,
                       max_clusters: int = 4) -> None:
    plot_xlim = DEFAULT_X_RANGE
    plot_ylim = DEFAULT_Y_RANGE
    plot_zlim = DEFAULT_Z_RANGE
    density = grid.density
    xy = np.sum(density, axis=2).T
    xz = np.sum(density, axis=1).T
    yz = np.sum(density, axis=0).T

    # 单独窗口：XY
    fig_xy, ax_xy = plt.subplots(figsize=(5, 4))
    im0 = ax_xy.imshow(xy, origin="lower", aspect="equal",
                       extent=[grid.x_range[0], grid.x_range[1], grid.y_range[0], grid.y_range[1]], cmap="jet")
    ax_xy.set_title("XY投影", fontsize=14)
    ax_xy.set_xlabel("X (mm)", fontsize=12)
    ax_xy.set_ylabel("Y (mm)", fontsize=12)
    ax_xy.tick_params(axis='both', labelsize=12)
    ax_xy.set_xlim(plot_xlim)
    ax_xy.set_ylim(plot_ylim)
    fig_xy.colorbar(im0, ax=ax_xy, shrink=0.8)
    fig_xy.tight_layout()
    out_xy = f"{prefix}_xy.png"
    fig_xy.savefig(out_xy, dpi=300)
    plt.close(fig_xy)
    print(f"保存 XY 投影到 {out_xy}")

    # 单独窗口：XZ
    fig_xz, ax_xz = plt.subplots(figsize=(5, 4))
    im1 = ax_xz.imshow(xz, origin="lower", aspect="equal",
                       extent=[grid.x_range[0], grid.x_range[1], grid.z_range[0], grid.z_range[1]], cmap="jet")
    ax_xz.set_title("XZ投影", fontsize=14)
    ax_xz.set_xlabel("X (mm)", fontsize=12)
    ax_xz.set_ylabel("Z (mm)", fontsize=12)
    ax_xz.tick_params(axis='both', labelsize=12)
    ax_xz.set_xlim(plot_xlim)
    ax_xz.set_ylim(plot_zlim)
    fig_xz.colorbar(im1, ax=ax_xz, shrink=0.8)
    fig_xz.tight_layout()
    out_xz = f"{prefix}_xz.png"
    fig_xz.savefig(out_xz, dpi=300)
    plt.close(fig_xz)
    print(f"保存 XZ 投影到 {out_xz}")

    # 单独窗口：YZ
    fig_yz, ax_yz = plt.subplots(figsize=(5, 4))
    im2 = ax_yz.imshow(yz, origin="lower", aspect="equal",
                       extent=[grid.y_range[0], grid.y_range[1], grid.z_range[0], grid.z_range[1]], cmap="jet")
    ax_yz.set_title("YZ投影", fontsize=14)
    ax_yz.set_xlabel("Y (mm)", fontsize=12)
    ax_yz.set_ylabel("Z (mm)", fontsize=12)
    ax_yz.tick_params(axis='both', labelsize=12)
    ax_yz.set_xlim(plot_ylim)
    ax_yz.set_ylim(plot_zlim)
    fig_yz.colorbar(im2, ax=ax_yz, shrink=0.8)
    fig_yz.tight_layout()
    out_yz = f"{prefix}_yz.png"
    fig_yz.savefig(out_yz, dpi=300)
    plt.close(fig_yz)
    print(f"保存 YZ 投影到 {out_yz}")


    # 顶部体素的 3D 视图（使用DBSCAN聚类去除离散点）
    # 先筛选出所有非零密度的体素
    nonzero_mask = density > 1e-10
    idxs = np.argwhere(nonzero_mask)
    
    if idxs.size > 0:
        centers = grid.origin + (idxs + 0.5) * grid.voxel_size
        values = density[nonzero_mask]
        
        # 使用DBSCAN聚类，eps为体素间的最大距离，min_samples为最小簇大小
        # 考虑密度作为权重：只对密度较高的点进行聚类
        density_threshold = np.percentile(values, density_percentile)
        valid_mask = values > density_threshold
        
        if np.sum(valid_mask) > 10:  # 至少有10个点
            centers_valid = centers[valid_mask]
            values_valid = values[valid_mask]
            
            # DBSCAN聚类：eps设为体素大小乘以因子
            dbscan_eps = grid.voxel_size * dbscan_eps_factor
            clustering = DBSCAN(eps=dbscan_eps, min_samples=dbscan_min_samples).fit(centers_valid)
            labels = clustering.labels_
            
            # 保留最大簇（标签不为-1的点，-1表示噪声）
            # 找出最大的簇
            unique_labels = set(labels)
            if -1 in unique_labels:
                unique_labels.remove(-1)  # 移除噪声标签
            
            if len(unique_labels) > 0:
                # 找出点数最多的几个簇
                label_counts = [(label, np.sum(labels == label)) for label in unique_labels]
                label_counts.sort(key=lambda x: x[1], reverse=True)
                
                # 保留最大的N个簇
                selected_labels = []
                for label, count in label_counts[:max_clusters]:
                    selected_labels.append(label)
                
                # 筛选出选中簇的点
                cluster_mask = np.isin(labels, selected_labels)
                centers_filtered = centers_valid[cluster_mask]
                values_filtered = values_valid[cluster_mask]
                
                # 如果点太多，随机采样
                if len(values_filtered) > 120000:
                    sel = np.random.choice(len(values_filtered), 120000, replace=False)
                    centers_filtered = centers_filtered[sel]
                    values_filtered = values_filtered[sel]
                
                # 绘制3D散点图
                fig3d = plt.figure(figsize=(5, 4))
                ax3d = fig3d.add_subplot(111, projection="3d")
                p = ax3d.scatter(centers_filtered[:, 0], centers_filtered[:, 1], centers_filtered[:, 2], 
                               c=values_filtered, cmap="turbo", s=1, alpha=0.8, edgecolors='none')
                ax3d.set_xlabel("X (mm)", fontsize=12)
                ax3d.set_ylabel("Y (mm)", fontsize=12)
                ax3d.set_zlabel("Z (mm)", fontsize=12)
                ax3d.set_xlim(plot_xlim)
                ax3d.set_ylim(plot_ylim)
                ax3d.set_zlim(plot_zlim)
                ax3d.set_box_aspect([1,1,1])
                ax3d.view_init(elev=50, azim=45)
                ax3d.set_title(f"密度体素图", fontsize=14, position=(0.6, 1.05))
                ax3d.tick_params(axis='both', labelsize=12)
                cbar = fig3d.colorbar(p, ax=ax3d, shrink=0.65)
                cbar.set_label('密度', fontsize=12)
                cbar.ax.tick_params(labelsize=12)
                
                # fig3d.tight_layout()
                out3d = f"{prefix}_3d.png"
                fig3d.savefig(out3d, dpi=500)
                plt.close(fig3d)
                print(f"保存3D视图到 {out3d} (DBSCAN聚类后保留{len(selected_labels)}个簇)")
            else:
                print("警告：聚类后无有效簇，跳过3D图")
        else:
            print("警告：有效点太少，跳过3D图")

    # 指标图拆分
    flat = density.ravel()
    fig_hist, ax_hist = plt.subplots(figsize=(5, 3))
    ax_hist.hist(flat, bins=100, log=True, color="steelblue", edgecolor="black")
    ax_hist.set_xlabel("密度", fontsize=12)
    ax_hist.set_ylabel("计数 (log)", fontsize=12)
    ax_hist.set_title("密度直方图", fontsize=14)
    ax_hist.tick_params(axis='both', labelsize=12)
    fig_hist.tight_layout()
    out_hist = f"{prefix}_density_hist.png"
    fig_hist.savefig(out_hist, dpi=300)
    plt.close(fig_hist)
    print(f"保存密度直方图到 {out_hist}")

    # 将目标 Z 转成索引
    z_indices = []
    for z in target_z_vals:
        # 按比例映射到索引
        ratio = (z - grid.z_range[0]) / (grid.z_range[1] - grid.z_range[0])
        idx = int(np.clip(round(ratio * (grid.nz - 1)), 0, grid.nz - 1))
        z_indices.append(idx)

    figz, axz = plt.subplots(2, 3, figsize=(12, 7))
    axz = axz.ravel()
    for i, (z_idx, z_val) in enumerate(zip(z_indices, target_z_vals)):
        slice_xy_full = density[:, :, z_idx].T
        slice_xy_i = slice_xy_full[::2, ::2]  # downsample to reduce pixel density
        vmin_i = np.percentile(slice_xy_i, 2)
        vmax_i = np.percentile(slice_xy_i, 98)
        im = axz[i].imshow(slice_xy_i, origin="lower", aspect="equal", cmap="jet",
                           extent=[grid.x_range[0], grid.x_range[1], grid.y_range[0], grid.y_range[1]],
                           vmin=0, vmax=50)
        axz[i].set_title(f"Z={z_val} mm", fontsize=14)
        axz[i].set_xlabel("X (mm)", fontsize=12)
        axz[i].set_ylabel("Y (mm)", fontsize=12)
        axz[i].tick_params(axis='both', labelsize=12)
        cbar = figz.colorbar(im, ax=axz[i], shrink=0.7)
        cbar.ax.tick_params(labelsize=12)

    # 若少于6个子图则隐藏空轴
    for j in range(len(target_z_vals), len(axz)):
        axz[j].axis('off')

    figz.suptitle("Z方向切片", fontsize=14)
    figz.tight_layout()
    outz = f"{prefix}_z_slices.png"
    figz.savefig(outz, dpi=500)
    plt.close(figz)
    print(f"保存 Z 切片到 {outz}")


def save_event_stats(events: List[EventData], rad_len_mm: float, prefix: str = "./picture/mlem") -> None:
    if not events:
        print("无事件可用于统计图。")
        return

    angles = np.array([e.angle_mrad for e in events], dtype=float)
    paths = np.array([e.total_length for e in events], dtype=float)
    
    x_over_x0_initial = paths / max(rad_len_mm, 1e-9)
    
    def angle_to_effective_x0(angle_mrad: float, momentum_mev: float) -> float:
        """
        根据散射角反推等效x/X0
        使用高地公式: theta0 = (13.6/pc) * sqrt(x/X0) * (1 + 0.038*ln(x/X0))
        简化为: x/X0 ≈ (theta0 * pc / 13.6)^2
        """
        theta0_rad = angle_mrad / 1000.0
        x_over_x0 = (theta0_rad * momentum_mev / HIGHLAND_K) ** 2
        return x_over_x0
    
    x_over_x0_effective = np.array([angle_to_effective_x0(a, DEFAULT_MUON_MOMENTUM) for a in angles])

    p90_angle = np.percentile(angles, 90)
    p50_angle = np.percentile(angles, 50)
    p90_x0_eff = np.percentile(x_over_x0_effective, 90)
    p50_x0_eff = np.percentile(x_over_x0_effective, 50)
    print(f"散射角: P50={p50_angle:.2f} mrad, P90={p90_angle:.2f} mrad")
    print(f"等效x/X0: P50={p50_x0_eff:.3f}, P90={p90_x0_eff:.3f}")

    fig_a, ax_a = plt.subplots(figsize=(5, 3))
    ax_a.hist(angles, bins=100, log=True, color="steelblue", edgecolor="black")
    ax_a.set_xlabel("散射角 (mrad)", fontsize=12)
    ax_a.set_ylabel("计数 (log)", fontsize=12)
    ax_a.set_title("散射角分布", fontsize=14)
    ax_a.tick_params(axis='both', labelsize=12)
    fig_a.tight_layout()
    out_a = f"{prefix}_angle_hist.png"
    fig_a.savefig(out_a, dpi=300)
    plt.close(fig_a)
    print(f"保存散射角直方图到 {out_a}")

    fig_x0, ax_x0 = plt.subplots(figsize=(5, 3))
    ax_x0.hist(x_over_x0_effective, bins=100, log=True, color="darkorange", edgecolor="black")
    ax_x0.set_xlabel("等效 $x/X_0$ ", fontsize=12)
    ax_x0.set_ylabel("计数 (log)", fontsize=12)
    ax_x0.set_title("等效 $x/X_0$ 分布", fontsize=14)
    ax_x0.tick_params(axis='both', labelsize=12)
    fig_x0.tight_layout()
    out_x0f = f"{prefix}_xOverX0_hist.png"
    fig_x0.savefig(out_x0f, dpi=300)
    plt.close(fig_x0)
    print(f"保存 x/X0 直方图到 {out_x0f}")

def build_events(raw_events: Dict[int, Dict[int, np.ndarray]], grid: VoxelGrid,
                 angle_min: float, angle_max: float) -> List[EventData]:
    events: List[EventData] = []
    for eid in sorted(raw_events):
        hits = raw_events[eid]
        if len(hits) < 4:
            continue
        p1, p2, p3, p4 = hits[1], hits[2], hits[3], hits[4]
        dir_in = p2 - p1
        dir_out_down = p4 - p3
        angle = scattering_angle_mrad(dir_in, dir_out_down)
        if angle < angle_min or angle > angle_max:
            continue

        # 传入线路：锚定在 L2，方向不变（向下）
        inc_contrib = trace_line(p2, dir_in, grid)
        # 出线：锚定在 L3，方向反转（向上）以穿过 ROI
        out_contrib = trace_line(p3, p3 - p4, grid)

        contributions: Dict[int, float] = defaultdict(float)
        for k, v in inc_contrib.items():
            contributions[k] += v
        for k, v in out_contrib.items():
            contributions[k] += v

        total_len = float(sum(contributions.values()))
        if total_len <= 1e-6:
            continue

        events.append(EventData(dict(contributions), angle, total_len))
    return events


def main() -> None:
    parser = argparse.ArgumentParser(description="EM-ML缪子成像，处理 build/output.csv 文件")
    parser.add_argument("--input", type=Path, default=Path("build/output.csv"), help="output.csv 文件路径")
    parser.add_argument("--voxel", type=float, default=10.0, help="体素大小 (mm)")
    parser.add_argument("--iterations", type=int, default=10, help="EM迭代次数")
    parser.add_argument("--max-events", type=int, default=50000, help="最大加载事件数 (近似; 需要完整事件)")
    parser.add_argument("--angle-min", type=float, default=1.0, help="最小散射角 (mrad)")
    parser.add_argument("--angle-max", type=float, default=200.0, help="最大散射角 (mrad)")
    parser.add_argument("--momentum", type=float, default=DEFAULT_MUON_MOMENTUM, help="Muon 动量 (MeV/c)")
    parser.add_argument("--radlen", type=float, default=DEFAULT_RADIATION_LENGTH, help="辐射长度 (mm)")
    parser.add_argument("--x-range", type=float, nargs=2, default=DEFAULT_X_RANGE, help="X 范围 (mm)")
    parser.add_argument("--y-range", type=float, nargs=2, default=DEFAULT_Y_RANGE, help="Y 范围 (mm)")
    parser.add_argument("--z-range", type=float, nargs=2, default=DEFAULT_Z_RANGE, help="Z 范围 (mm)")
    parser.add_argument("--density-percentile", type=float, default=90.0,
                        help="DBSCAN前筛选的密度分位数 (默认90)")
    parser.add_argument("--dbscan-eps-factor", type=float, default=1.0,
                        help="DBSCAN的eps相对于体素大小的倍数 (默认1.0)")
    parser.add_argument("--dbscan-min-samples", type=int, default=3,
                        help="DBSCAN的最小样本数 (默认3)")
    parser.add_argument("--max-clusters", type=int, default=4,
                        help="最多保留的聚类簇数 (默认4)")
    parser.add_argument("--no-plot", action="store_true", help="跳过保存投影图")
    args = parser.parse_args()

    grid = VoxelGrid(tuple(args.x_range), tuple(args.y_range), tuple(args.z_range), args.voxel)
    print(f"网格尺寸: {grid.nx} x {grid.ny} x {grid.nz} = {grid.total_voxels()} 个体素")

    raw_events = parse_output(args.input, max_events=args.max_events)
    print(f"从 {args.input} 加载了 {len(raw_events)} 个事件")

    events = build_events(raw_events, grid, args.angle_min, args.angle_max)
    print(f"使用 {len(events)} 个事件经过角度/几何裁剪")

    recon = EMLReconstructor(grid, events, args.momentum, args.radlen)
    recon.run(args.iterations)
    
    if not args.no_plot:
        save_density_plots(grid, prefix="./picture/mlem",
                           density_percentile=args.density_percentile,
                           dbscan_eps_factor=args.dbscan_eps_factor,
                           dbscan_min_samples=args.dbscan_min_samples,
                           max_clusters=args.max_clusters)
        save_event_stats(events, args.radlen)

if __name__ == "__main__":
    main()
