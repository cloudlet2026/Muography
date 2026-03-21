'''https://github.com/BehnamManafi/Poca_algorithm'''
import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib import font_manager
from mpl_toolkits.mplot3d import Axes3D

# ---------------- 可调过滤参数 ----------------
# 最小散射角阈值（毫弧度），用于剔除小角度散射事件
SCATTER_MIN_MRAD = 1.0
# 是否应用PoCA点Z窗口过滤（将PoCA限定在物体厚度范围附近）
USE_Z_WINDOW = False
# 物体厚度为5cm -> 以中心为0，±120mm作为默认窗口
POCA_Z_WINDOW_MM = 120.0
# ------------------------------------------------
    # 设置默认字体
font_manager.fontManager.addfont("fonts/TimesSimSunRegular.ttf")
plt.rcParams["font.family"] = "TimesSimSun" 

# ---------------- 几何参数 ----------------
# 探测器层z坐标（毫米）
Z1 = 800.0   # 模块1
Z2 = 600.0    # 模块2
Z3 = -600.0   # 模块3
Z4 = -800.0  # 模块4
# ------------------------------------------

def parse_output_file(filename):
    """
    解析 output.csv 文件
    返回字典: {event_id: {layer_id: coordinates}}
    """
    events = {}
    with open(filename, 'r') as f:
        header = f.readline().strip()
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split(',')
            if len(parts) >= 4:
                event_id = int(parts[0])
                layer_id = int(parts[1])
                x = float(parts[2])
                y = float(parts[3])
                
                # 根据moduleId设置z坐标
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
                events[event_id][layer_id] = np.array([x, y, z])
    
    return events

def poca_algorithm(line1_start, line1_dir, line2_start, line2_dir, clamp_to_segment=True):
    """
    计算两条异面直线的最接近点
    直线1: P1(t) = line1_start + t * line1_dir
    直线2: P2(s) = line2_start + s * line2_dir
    
    参数:
        clamp_to_segment: 是否将参数限制在[0,1]范围（线段而非无限直线）
    
    返回: (最接近点, 距离, 上层线上的点, 下层线上的点, t, s)
    """
    # 保留原始向量长度
    d1 = line1_dir
    d2 = line2_dir
    
    # 两条直线的起点连接向量
    w = line1_start - line2_start
    
    # 计算参数
    a = np.dot(d1, d1)
    b = np.dot(d1, d2)
    c = np.dot(d2, d2)
    d = np.dot(d1, w)
    e = np.dot(d2, w)
    
    # 分母
    denom = a * c - b * b
    
    if abs(denom) < 1e-10:
        # 两直线平行
        t = 0.5
        s = 0.5
    else:
        # 最接近点的参数
        t = (b * e - c * d) / denom
        s = (a * e - b * d) / denom
    
    # 如果启用约束，将参数限制在[0,1]范围内
    if clamp_to_segment:
        t = np.clip(t, 0, 1)
        s = np.clip(s, 0, 1)
    
    # 两条直线上的最接近点
    poca_line1 = line1_start + t * d1
    poca_line2 = line2_start + s * d2
    
    # 计算最小距离
    distance = np.linalg.norm(poca_line1 - poca_line2)
    
    # 计算最接近点的坐标（两点的中点）
    closest_approach = (poca_line1 + poca_line2) / 2.0
    
    return closest_approach, distance, poca_line1, poca_line2, t, s

def check_poca_validity(poca_point, z_min=-120.0, z_max=120.0):
    """
    检查 PoCA 点是否有效
    返回: (是否有效, 标记字符串)
    """
    z_coord = poca_point[2]
    if z_min <= z_coord <= z_max:
        return True
    else:
        return False

def filter_isolated_points(points, radius=10.0, min_neighbors=10):
    """
    使用三维网格的邻域统计过滤离散点：
    - 将空间划分为边长为 `radius` 的立方体网格
    - 仅在本格子及其26个相邻格子中统计邻居，避免 O(N^2) 内存/时间

    参数:
        points: (N, 3) 点集
        radius: 邻域半径（mm）
        min_neighbors: 至少需要的邻居数量（不含自身）

    返回:
        (filtered_points, mask)
    """
    if points is None or len(points) == 0:
        return points, np.array([], dtype=bool)

    r2 = radius * radius
    # 计算每个点所在的网格坐标
    cells = np.floor(points / radius).astype(int)
    cell_dict = {}
    for idx, cell in enumerate(map(tuple, cells)):
        if cell not in cell_dict:
            cell_dict[cell] = []
        cell_dict[cell].append(idx)

    # 预生成邻居偏移（3x3x3 共27个）
    offsets = [(dx, dy, dz) for dx in (-1, 0, 1)
                          for dy in (-1, 0, 1)
                          for dz in (-1, 0, 1)]

    keep_mask = np.zeros(len(points), dtype=bool)

    for i, cell in enumerate(map(tuple, cells)):
        cx, cy, cz = cell
        # 收集邻近格子中的候选点索引
        candidates = []
        for dx, dy, dz in offsets:
            nbr_cell = (cx + dx, cy + dy, cz + dz)
            if nbr_cell in cell_dict:
                candidates.extend(cell_dict[nbr_cell])

        if not candidates:
            continue

        cand_idx = np.array(candidates, dtype=int)
        cand_pts = points[cand_idx]
        diff = cand_pts - points[i]
        d2 = np.einsum('ij,ij->i', diff, diff)
        # 邻居数（不计自身：距离为0）
        neighbor_count = int(np.sum(d2 <= r2)) - 1
        if neighbor_count >= min_neighbors:
            keep_mask[i] = True

    return points[keep_mask], keep_mask

def scattering_angle_mrad(dir_top, dir_bottom_downward):
    """
    计算散射角（毫弧度）。
    参数:
        dir_top: 上层段方向（朝向下方）
        dir_bottom_downward: 下层段方向（同样朝向下方）
    返回:
        angle_mrad: 散射角（mrad）
    """
    u1 = dir_top / (np.linalg.norm(dir_top) + 1e-12)
    u2 = dir_bottom_downward / (np.linalg.norm(dir_bottom_downward) + 1e-12)
    cos_theta = np.clip(np.dot(u1, u2), -1.0, 1.0)
    theta = np.arccos(cos_theta)
    return float(theta * 1000.0)

# 解析 output.csv 文件
events = parse_output_file('./build/output.csv')

print("=" * 70)
print("PoCA 算法 - 事件分析")
print("=" * 70)

# 统计信息
poca_stats = {'有效': 0, '无效': 0}
filtered_small_angle = 0
filtered_z_window = 0

# 存储所有PoCA点用于绘图
poca_points = []
scatter_angles = []  # 散射角（mrad）
poca_distances = []  # PoCA距离（mm）

# 处理每个事件
for event_id in sorted(events.keys()):
    event_data = events[event_id]
    
    # 检查是否有所有4个层的数据
    if len(event_data) < 4:
        print(f"\n事件 {event_id}: 数据不完整 (只有 {len(event_data)} 层)")
        continue
    
    # 提取各层坐标
    l1_coord = event_data[1]  # L1 (顶部, z=300mm)
    l2_coord = event_data[2]  # L2 (z=300mm)
    l3_coord = event_data[3]  # L3 (z=-300mm)
    l4_coord = event_data[4]  # L4 (底部, z=-300mm)
    
    # 计算上层线 (L1 -> L2)
    line1_start = l1_coord
    line1_end = l2_coord
    line1_dir = line1_end - line1_start
    
    # 计算下层线 (L4 -> L3)
    line2_start = l4_coord
    line2_end = l3_coord
    line2_dir = line2_end - line2_start
    
    # 计算 PoCA
    closest_approach, distance, poca_line1, poca_line2, t, s = poca_algorithm(
        line1_start, line1_dir, line2_start, line2_dir, clamp_to_segment=False
    )
    
    # 计算散射角：
    # 顶部段采用 L1->L2 （向下），底部段应与其同向，故使用 L4->L3 的反向即 L3->L4
    top_dir_down = line1_dir
    bottom_dir_down = l4_coord - l3_coord  # 与下层真实通过方向一致（向下）
    angle_mrad = scattering_angle_mrad(top_dir_down, bottom_dir_down)

    # 过滤：散射角过小
    if angle_mrad < SCATTER_MIN_MRAD:
        filtered_small_angle += 1
        continue

    # 过滤：Z窗口（可选）
    if USE_Z_WINDOW and (abs(closest_approach[2]) > POCA_Z_WINDOW_MM):
        filtered_z_window += 1
        continue

    # 检查 PoCA 点的有效性（基础范围）
    is_valid = check_poca_validity(closest_approach)
    
    # 存储PoCA点及指标
    poca_points.append(closest_approach.copy())
    scatter_angles.append(angle_mrad)
    poca_distances.append(distance)
    
    # 检查参数是否在合理范围内
    t_valid = 0 <= t <= 1
    s_valid = 0 <= s <= 1
    param_warning = "" if (t_valid and s_valid) else " [参数超限]"
    
    # 确定输出标记
    if is_valid:
        marker = "[有效]"
        poca_stats['有效'] += 1
    else:
        marker = "[无效]"
        poca_stats['无效'] += 1
    
    # 只显示前10个事件
    if event_id <= 10:
        print(f"\n事件 {event_id}{param_warning}:")
        print(f"  PoCA 点: {closest_approach}")
        print(f"  PoCA z: {closest_approach[2]:.2f} mm  | 散射角: {angle_mrad:.2f} mrad")


print("\n" + "=" * 70)
print("统计信息:")
print(f"  有效 PoCA 点: {poca_stats['有效']} 个")
print(f"  无效 PoCA 点: {poca_stats['无效']} 个")
print(f"  过滤(角度< {SCATTER_MIN_MRAD} mrad): {filtered_small_angle} 个")
if USE_Z_WINDOW:
    print(f"  过滤(|z|> {POCA_Z_WINDOW_MM} mm): {filtered_z_window} 个")
print(f"  总计: {poca_stats['有效'] + poca_stats['无效']} 个事件")
print("=" * 70)

# 绘图
plot_xlim = (-600,600)
plot_ylim = (-600,600)
plot_zlim = (-600,600)
if len(poca_points) > 0:
    poca_array = np.array(poca_points)

    # 邻域密度过滤，剔除比较离散的点
    min_neighbors=8
    radius=40.0
    filtered_points, _mask = filter_isolated_points(poca_array, radius=radius, min_neighbors=min_neighbors)
    if filtered_points is None or len(filtered_points) == 0:
        plot_points = poca_array
        print("密度过滤后无点保留，使用原始点绘图。")
    else:
        plot_points = filtered_points
        removed = len(poca_array) - len(plot_points)
        print(f"密度过滤：移除 {removed} 个离散点（半径{radius}mm, 至少 {min_neighbors} 邻居）。")
        
    range2d = (plot_xlim[0], plot_xlim[1])
    range3d = [[plot_xlim[0], plot_xlim[1]], [plot_ylim[0], plot_ylim[1]]]
    # 图1: 三维散点图
    fig1 = plt.figure(figsize=(5, 4))
    ax1 = fig1.add_subplot(111, projection='3d')
    
    # 使用单一、平静的颜色来避免繁杂的渐变
    scatter_color = 'green'
    scatter = ax1.scatter(plot_points[:, 0], plot_points[:, 1], plot_points[:, 2],
                         color=scatter_color, s=1 , alpha=1.0)
    
    # 绘制散射体位置平面
    xx, yy = np.meshgrid([-400, 400], [-400, 400])
    ax1.view_init(elev=25, azim=60)
    ax1.set_xlabel('X (mm)',fontsize=12)
    ax1.set_ylabel('Y (mm)',fontsize=12)
    ax1.set_zlabel('Z (mm)',fontsize=12)
    ax1.set_title('PoCA点三维分布',fontsize=14,position=(0.6, 1.05))
    ax1.set_xlim([plot_xlim[0], plot_xlim[1]])
    ax1.set_ylim([plot_ylim[0], plot_ylim[1]])
    ax1.set_zlim([plot_zlim[0], plot_zlim[1]])
    ax1.set_box_aspect([1,1,1])
    ax1.view_init(elev=25, azim=45)
    ax1.tick_params(axis='both', which='major', labelsize=12)
    # plt.tight_layout()
    plt.savefig('./picture/poca_3d_visualization.png', dpi=300)
    print("\n3D图已保存到: poca_3d_visualization.png")
    
    # 图2: 二维投影图（XY平面）
    fig2 = plt.figure(figsize=(5, 4))
    ax2 = fig2.add_subplot(111)
    scatter2 = ax2.scatter(plot_points[:, 0], plot_points[:, 1],
                          color=scatter_color, s=1, alpha=1.0)
    
    # 绘制散射体边界（在XY平面的投影）
    # square = plt.Rectangle((-50, -50), 100, 100, fill=False, color='red', linewidth=2, 
    #                      linestyle='--', label='散射体区域 (±5cm)')
    # ax2.add_patch(square)
    
    ax2.set_xlabel('X (mm)',fontsize=12)
    ax2.set_ylabel('Y (mm)',fontsize=12)
    ax2.set_title('PoCA点二维分布 (XY投影)',fontsize=14)
    ax2.set_xlim([plot_xlim[0], plot_xlim[1]])
    ax2.set_ylim([plot_ylim[0], plot_ylim[1]])
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    # ax2.legend(fontsize=12)
    
    plt.tight_layout()
    plt.savefig('./picture/poca_2d_projection.png', dpi=300)
    print("2D投影图已保存到: poca_2d_projection.png")
    
    # ========== 创建质量指标图 ==========
    print("\n正在生成质量指标图...")
    
    # 图3: 密度热图（2D直方图）
    fig3 = plt.figure(figsize=(5, 4))
    ax = fig3.add_subplot(111)
    h, xedges, yedges = np.histogram2d(plot_points[:, 0], plot_points[:, 1], bins=80, range=range3d)
    im = ax.imshow(h.T, origin='lower', aspect='equal',
                   extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                   cmap='jet')
    ax.set_xlabel('X (mm)',fontsize=12)
    ax.set_ylabel('Y (mm)',fontsize=12)
    ax.set_title('密度热图',fontsize=14)
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('事件计数', fontsize=12)
    cbar.ax.tick_params(labelsize=12)
    # square = plt.Rectangle((-200, -200), 400, 400, fill=False, color='orange', linewidth=2, linestyle='--', label='散射体区域 (±20cm)')
    # ax.add_patch(square)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('./picture/poca_density_heatmap.png', dpi=300)
    print("密度热图已保存到: poca_density_heatmap.png")
    
    # 图4: X轴投影
    fig4 = plt.figure(figsize=(5, 3))
    ax = fig4.add_subplot(111)
    x_proj, x_bins = np.histogram(plot_points[:, 0], bins=80, range=range2d)
    x_centers = (x_bins[:-1] + x_bins[1:]) / 2
    ax.bar(x_centers, x_proj, width=x_bins[1]-x_bins[0], color='blue', alpha=0.7, edgecolor='black')
    ax.set_xlabel('X (mm)',fontsize=12)
    ax.set_ylabel('事件计数',fontsize=12)
    ax.set_title('X轴投影分布',fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('./picture/poca_x_projection.png', dpi=300)
    print("X轴投影已保存到: poca_x_projection.png")
    
    # 图5: Y轴投影
    fig5 = plt.figure(figsize=(5, 3))
    ax = fig5.add_subplot(111)
    y_proj, y_bins = np.histogram(plot_points[:, 1], bins=80, range=range2d)
    y_centers = (y_bins[:-1] + y_bins[1:]) / 2
    ax.bar(y_centers, y_proj, width=y_bins[1]-y_bins[0], color='red', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Y (mm)',fontsize=12)
    ax.set_ylabel('事件计数',fontsize=12)
    ax.set_title('Y轴投影分布',fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('./picture/poca_y_projection.png', dpi=300)
    print("Y轴投影已保存到: poca_y_projection.png")
    
    # 图6: Z坐标分布
    fig6 = plt.figure(figsize=(5, 3))
    ax = fig6.add_subplot(111)
    z_proj, z_bins = np.histogram(plot_points[:, 2], bins=80, range=plot_zlim)
    z_centers = (z_bins[:-1] + z_bins[1:]) / 2
    ax.bar(z_centers, z_proj, width=z_bins[1]-z_bins[0], color='green', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Z (mm)',fontsize=12)
    ax.set_ylabel('事件计数',fontsize=12)
    ax.set_title('Z坐标分布',fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('./picture/poca_z_distribution.png', dpi=300)
    print("Z坐标分布已保存到: poca_z_distribution.png")
    
    # 图7: 散射角分布
    fig7 = plt.figure(figsize=(5, 3))
    ax = fig7.add_subplot(111)
    # 只绘制通过过滤的点对应的散射角
    valid_angles = np.array(scatter_angles)
    ax.hist(valid_angles, bins=100, color='purple', alpha=0.7, edgecolor='black')
    ax.set_xlabel('散射角 (mrad)',fontsize=12)
    ax.set_ylabel('事件计数',fontsize=12)
    ax.set_title('散射角分布',fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('./picture/poca_scatter_angle.png', dpi=300)
    print("散射角分布已保存到: poca_scatter_angle.png")
    
    # 图8: PoCA距离分布
    fig8 = plt.figure(figsize=(5, 3))
    ax = fig8.add_subplot(111)
    valid_distances = np.array(poca_distances)
    ax.hist(valid_distances, bins=100, color='brown', alpha=0.7, edgecolor='black')
    ax.set_xlabel('PoCA距离 (mm)',fontsize=12)
    ax.set_ylabel('事件计数',fontsize=12)
    ax.set_title('PoCA最小距离分布',fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('./picture/poca_distance.png', dpi=300)
    print("PoCA距离分布已保存到: poca_distance.png")
    
    # 显示所有图形
    # plt.show()
    
    # 打印详细统计
    print("\n========== 详细统计信息 ==========")
    print(f"过滤后有效点数: {len(plot_points)}")
    print(f"\nPoCA点空间分布:")
    print(f"  X范围: [{plot_points[:, 0].min():.2f}, {plot_points[:, 0].max():.2f}] mm")
    print(f"  Y范围: [{plot_points[:, 1].min():.2f}, {plot_points[:, 1].max():.2f}] mm")
    print(f"  Z范围: [{plot_points[:, 2].min():.2f}, {plot_points[:, 2].max():.2f}] mm")
    print(f"  X中心: {plot_points[:, 0].mean():.2f} ± {plot_points[:, 0].std():.2f} mm")
    print(f"  Y中心: {plot_points[:, 1].mean():.2f} ± {plot_points[:, 1].std():.2f} mm")
    print(f"  Z中心: {plot_points[:, 2].mean():.2f} ± {plot_points[:, 2].std():.2f} mm")
    print(f"\n散射角统计:")
    print(f"  范围: [{valid_angles.min():.2f}, {valid_angles.max():.2f}] mrad")
    print(f"  均值: {valid_angles.mean():.2f} ± {valid_angles.std():.2f} mrad")
    print(f"  中位数: {np.median(valid_angles):.2f} mrad")
    print(f"\nPoCA距离统计:")
    print(f"  范围: [{valid_distances.min():.2f}, {valid_distances.max():.2f}] mm")
    print(f"  均值: {valid_distances.mean():.2f} ± {valid_distances.std():.2f} mm")
    print(f"  中位数: {np.median(valid_distances):.2f} mm")
    print("=" * 70)

