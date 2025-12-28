"""
重新生成粒子追踪可视化图，添加底图
"""

import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import pickle
from scipy.ndimage import gaussian_filter
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# 文件路径
TRAJECTORY_FILE = r'E:\HK_model\202206_runtime\particle_tracking\optimized_results\trajectory.pkl'
SHP_FILE = r'E:\HK_model\202206_runtime\particle_tracking\PRD.shp'
OUTPUT_DIR = r'E:\HK_model\202206_runtime\particle_tracking\optimized_results'

print("=" * 60)
print("重新生成粒子追踪可视化图")
print("=" * 60)

# 加载轨迹数据
print("\n[1/3] 加载轨迹数据...")
with open(TRAJECTORY_FILE, 'rb') as f:
    data = pickle.load(f)

trajectory = data['trajectory']
source_ids = data['source_ids']
sources = data['sources']
n_particles = trajectory.shape[1]

print(f"[OK] 轨迹形状: {trajectory.shape}")
print(f"[OK] 粒子数: {n_particles}")
print(f"[OK] 时间点数: {trajectory.shape[0]}")

# 加载底图
print("\n[2/3] 加载底图...")
try:
    gdf = gpd.read_file(SHP_FILE)
    bounds = gdf.total_bounds  # [minx, miny, maxx, maxy]
    print(f"[OK] 底图范围: {bounds}")
    print(f"    经度: {bounds[0]:.4f} - {bounds[2]:.4f}")
    print(f"    纬度: {bounds[1]:.4f} - {bounds[3]:.4f}")
except Exception as e:
    print(f"[ERROR] 加载底图失败: {e}")
    gdf = None
    # 使用默认范围
    bounds = [113.0, 21.5, 114.5, 23.0]

# ==================== 重新生成密度图 ====================
print("\n[3/3] 生成密度图...")

# 获取最终位置
final_pos = trajectory[-1]
valid_mask = ~(np.isnan(final_pos[:, 0]) | np.isnan(final_pos[:, 1]))
final_lons = final_pos[valid_mask, 0]
final_lats = final_pos[valid_mask, 1]

# 限制在底图范围内
in_bounds = ((final_lons >= bounds[0]) & (final_lons <= bounds[2]) &
             (final_lats >= bounds[1]) & (final_lats <= bounds[3]))
final_lons = final_lons[in_bounds]
final_lats = final_lats[in_bounds]

print(f"底图范围内的粒子数: {len(final_lons)}")

# 创建密度图
fig, ax = plt.subplots(figsize=(16, 14))

# 绘制底图
if gdf is not None:
    gdf.plot(ax=ax, color='lightgray', edgecolor='black', linewidth=0.5, alpha=0.5)

# 2D直方图 - 使用底图范围
lon_bins = np.linspace(bounds[0], bounds[2], 200)
lat_bins = np.linspace(bounds[1], bounds[3], 200)

H, xedges, yedges = np.histogram2d(final_lons, final_lats, bins=[lon_bins, lat_bins])
H_smooth = gaussian_filter(H.T, sigma=2)

# 绘制密度
im = ax.imshow(
    H_smooth,
    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    origin='lower',
    cmap='hot',
    aspect='auto',
    alpha=0.7
)

# 标记排放口
ax.scatter(
    sources['longitude'],
    sources['latitude'],
    c='cyan',
    s=300,
    marker='*',
    edgecolors='white',
    linewidths=2,
    label='Emission Sources',
    zorder=10
)

for idx, source in sources.iterrows():
    ax.annotate(
        source['Point'],
        (source['longitude'], source['latitude']),
        xytext=(8, 8),
        textcoords='offset points',
        fontsize=11,
        color='white',
        weight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='blue', alpha=0.8),
        zorder=11
    )

plt.colorbar(im, ax=ax, label='Particle Density', fraction=0.046, pad=0.04)
ax.set_xlabel('Longitude (°E)', fontsize=13)
ax.set_ylabel('Latitude (°N)', fontsize=13)
ax.set_title(
    'Particle Density Distribution\n2022-06-01 to 2022-07-11',
    fontsize=15,
    weight='bold'
)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xlim(bounds[0], bounds[2])
ax.set_ylim(bounds[1], bounds[3])

density_file = f'{OUTPUT_DIR}/particle_density_with_basemap.png'
plt.savefig(density_file, dpi=300, bbox_inches='tight')
print(f"[OK] 密度图已保存: {density_file}")
plt.close()

# ==================== 重新生成轨迹图 ====================
print("\n生成轨迹图...")

fig, ax = plt.subplots(figsize=(16, 14))

# 绘制底图
if gdf is not None:
    gdf.plot(ax=ax, color='lightgray', edgecolor='black', linewidth=0.5, alpha=0.5)

# 绘制轨迹（采样）
n_plot = min(100, n_particles)
sample_indices = np.random.choice(n_particles, n_plot, replace=False)

for idx in sample_indices:
    traj = trajectory[:, idx, :]
    valid = ~(np.isnan(traj[:, 0]) | np.isnan(traj[:, 1]))

    # 限制在底图范围内
    traj_in_bounds = traj[valid]
    in_range = ((traj_in_bounds[:, 0] >= bounds[0]) & (traj_in_bounds[:, 0] <= bounds[2]) &
                (traj_in_bounds[:, 1] >= bounds[1]) & (traj_in_bounds[:, 1] <= bounds[3]))

    if np.any(in_range):
        ax.plot(traj_in_bounds[:, 0], traj_in_bounds[:, 1], alpha=0.4, linewidth=0.8)

# 标记排放口
ax.scatter(
    sources['longitude'],
    sources['latitude'],
    c='red',
    s=300,
    marker='*',
    edgecolors='white',
    linewidths=2,
    label='Emission Sources',
    zorder=10
)

for idx, source in sources.iterrows():
    ax.annotate(
        source['Point'],
        (source['longitude'], source['latitude']),
        xytext=(8, 8),
        textcoords='offset points',
        fontsize=11,
        color='black',
        weight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8),
        zorder=11
    )

ax.set_xlabel('Longitude (°E)', fontsize=13)
ax.set_ylabel('Latitude (°N)', fontsize=13)
ax.set_title(
    f'Particle Trajectories (Sample: {n_plot} particles)\n2022-06-01 to 2022-07-11',
    fontsize=15,
    weight='bold'
)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_xlim(bounds[0], bounds[2])
ax.set_ylim(bounds[1], bounds[3])

traj_file = f'{OUTPUT_DIR}/particle_trajectories_with_basemap.png'
plt.savefig(traj_file, dpi=300, bbox_inches='tight')
print(f"[OK] 轨迹图已保存: {traj_file}")
plt.close()

print("\n" + "=" * 60)
print("完成！")
print("=" * 60)
print(f"结果保存在: {OUTPUT_DIR}")
print(f"- 密度图: particle_density_with_basemap.png")
print(f"- 轨迹图: particle_trajectories_with_basemap.png")
print("=" * 60)
