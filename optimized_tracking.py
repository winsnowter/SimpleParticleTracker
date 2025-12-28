"""
优化版粒子追踪 - 专门处理大型 Delft3D FM 文件
优化策略：
1. 减少粒子数量
2. 增加时间步长和保存间隔
3. 分段运行避免内存溢出
4. 实时保存避免数据丢失
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from datetime import datetime, timedelta
import netCDF4 as nc
from tqdm import tqdm
import pickle
import os
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

# ==================== 配置参数 ====================
MAP_FILE = r'E:\HK_model\202206_runtime\202206_June_run.dsproj_data\HK-DFM\output\HK-DFM_map.nc'
SOURCE_FILE = r'E:\HK_model\202206_runtime\point.xlsx'
OUTPUT_DIR = r'E:\HK_model\202206_runtime\particle_tracking\optimized_results'

# 时间范围
START_DATE = datetime(2022, 6, 1)
END_DATE = datetime(2022, 7, 11)
BASE_DATE = datetime(2022, 1, 1)

# 优化参数
N_PARTICLES_PER_SOURCE = 100      # 减少粒子数（之前是1000）
DT = 3600                          # 时间步长1小时（之前30分钟）
SAVE_INTERVAL = 48                 # 每48步保存一次（2天保存一次）
DIFFUSION_COEF = 0.5              # 扩散系数
VELOCITY_LAYER = 'average'         # 深度平均流速
CHECKPOINT_INTERVAL = 200          # 每200步保存检查点

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("优化版粒子追踪系统")
print("=" * 60)

# ==================== 1. 加载数据 ====================
print("\n[1/4] 加载数据...")
print("正在打开 NetCDF 文件...")
ds = nc.Dataset(MAP_FILE, 'r')

# 读取网格
face_x = ds.variables['mesh2d_face_x'][:]
face_y = ds.variables['mesh2d_face_y'][:]
times = ds.variables['time'][:]

# 计算时间范围
start_seconds = (START_DATE - BASE_DATE).total_seconds()
end_seconds = (END_DATE - BASE_DATE).total_seconds()
start_idx = np.argmin(np.abs(times - start_seconds))
end_idx = np.argmin(np.abs(times - end_seconds))

print(f"[OK] 网格单元数: {len(face_x)}")
print(f"[OK] 总时间步数: {len(times)}")
print(f"[OK] 使用时间步范围: {start_idx} - {end_idx} ({end_idx-start_idx+1} 步)")
print(f"[OK] 实际开始时间: {BASE_DATE + timedelta(seconds=float(times[start_idx]))}")
print(f"[OK] 实际结束时间: {BASE_DATE + timedelta(seconds=float(times[end_idx]))}")

# 建立空间索引
print("正在建立空间索引...")
kdtree = cKDTree(np.column_stack([face_x, face_y]))
print("[OK] 空间索引建立完成")

# 加载排放口
sources = pd.read_excel(SOURCE_FILE)
print(f"[OK] 排放口数量: {len(sources)}")

# ==================== 2. 释放粒子 ====================
print("\n[2/4] 释放粒子...")
particles_list = []
source_ids_list = []

for idx, source in sources.iterrows():
    x0 = source['longitude']
    y0 = source['latitude']
    # 小范围随机分布（约±100m）
    x_p = x0 + np.random.normal(0, 0.001, N_PARTICLES_PER_SOURCE)
    y_p = y0 + np.random.normal(0, 0.001, N_PARTICLES_PER_SOURCE)
    particles_list.append(np.column_stack([x_p, y_p]))
    source_ids_list.extend([idx] * N_PARTICLES_PER_SOURCE)

particles = np.vstack(particles_list)
source_ids = np.array(source_ids_list)
n_particles = len(particles)

print(f"[OK] 总粒子数: {n_particles}")

# ==================== 3. 粒子追踪 ====================
print(f"\n[3/4] 开始粒子追踪...")
print(f"时间步长: {DT}秒 ({DT/3600}小时)")
print(f"保存间隔: 每{SAVE_INTERVAL}步")
print(f"检查点间隔: 每{CHECKPOINT_INTERVAL}步")

n_steps = end_idx - start_idx + 1
n_saved = n_steps // SAVE_INTERVAL + 1

# 初始化存储
trajectory = np.zeros((n_saved, n_particles, 2))
trajectory[0] = particles.copy()
time_indices = np.zeros(n_saved, dtype=int)
time_indices[0] = start_idx

# 当前位置
pos = particles.copy()
sigma = np.sqrt(2 * DIFFUSION_COEF * DT)
save_counter = 0

print(f"\n预计总步数: {n_steps}")
print(f"预计保存点数: {n_saved}")
print("\n开始追踪...")

for step in tqdm(range(n_steps), desc="追踪进度"):
    current_time_idx = start_idx + step

    if current_time_idx >= len(times):
        print(f"\n已达到时间序列末端")
        break

    # 读取流速
    try:
        u_all = ds.variables['mesh2d_ucx'][current_time_idx, :, :]
        v_all = ds.variables['mesh2d_ucy'][current_time_idx, :, :]
    except Exception as e:
        print(f"\n时间步{current_time_idx}读取失败: {e}")
        break

    # 为每个粒子计算流速
    velocities = np.zeros((n_particles, 2))
    for i in range(n_particles):
        try:
            dist, idx = kdtree.query([pos[i, 0], pos[i, 1]], k=1)
            if dist > 0.05:  # 距离太远
                continue

            # 深度平均流速
            if VELOCITY_LAYER == 'average':
                u_layers = u_all[idx, :]
                v_layers = v_all[idx, :]
                if np.ma.is_masked(u_layers):
                    u_valid = u_layers[~u_layers.mask]
                    v_valid = v_layers[~v_layers.mask]
                else:
                    u_valid = u_layers
                    v_valid = v_layers

                if len(u_valid) > 0 and len(v_valid) > 0:
                    velocities[i, 0] = u_valid.mean()
                    velocities[i, 1] = v_valid.mean()
        except:
            continue

    # 坐标转换
    lat_rad = np.radians(pos[:, 1])
    meters_per_deg_lon = 111320 * np.cos(lat_rad)
    meters_per_deg_lat = 111320

    # 平流项
    dlon_dt = velocities[:, 0] / meters_per_deg_lon
    dlat_dt = velocities[:, 1] / meters_per_deg_lat
    advection = np.column_stack([dlon_dt, dlat_dt]) * DT

    # 扩散项
    diffusion_m = np.random.normal(0, sigma, (n_particles, 2))
    diffusion_lon = diffusion_m[:, 0] / meters_per_deg_lon
    diffusion_lat = diffusion_m[:, 1] / meters_per_deg_lat
    diffusion = np.column_stack([diffusion_lon, diffusion_lat])

    # 更新位置
    pos += advection + diffusion

    # 保存轨迹
    if (step + 1) % SAVE_INTERVAL == 0:
        save_counter += 1
        if save_counter < n_saved:
            trajectory[save_counter] = pos.copy()
            time_indices[save_counter] = current_time_idx

    # 保存检查点
    if (step + 1) % CHECKPOINT_INTERVAL == 0:
        checkpoint_file = os.path.join(OUTPUT_DIR, f'checkpoint_step_{step+1}.pkl')
        with open(checkpoint_file, 'wb') as f:
            pickle.dump({
                'pos': pos.copy(),
                'step': step + 1,
                'save_counter': save_counter
            }, f)
        print(f"\n检查点已保存: step {step+1}")

# 裁剪到实际保存的步数
trajectory = trajectory[:save_counter+1]
time_indices = time_indices[:save_counter+1]

print(f"\n[OK] 追踪完成！实际保存了 {save_counter+1} 个时间点")

# 保存结果
results_file = os.path.join(OUTPUT_DIR, 'trajectory.pkl')
with open(results_file, 'wb') as f:
    pickle.dump({
        'trajectory': trajectory,
        'time_indices': time_indices,
        'source_ids': source_ids,
        'sources': sources,
        'time': times,
        'face_x': face_x,
        'face_y': face_y
    }, f)
print(f"[OK] 结果已保存: {results_file}")

# ==================== 4. 生成密度图 ====================
print("\n[4/4] 生成密度图...")

# 最终位置
final_pos = trajectory[-1]
valid_mask = ~(np.isnan(final_pos[:, 0]) | np.isnan(final_pos[:, 1]))
final_lons = final_pos[valid_mask, 0]
final_lats = final_pos[valid_mask, 1]

# 创建密度图
fig, ax = plt.subplots(figsize=(14, 12))

# 2D直方图
lon_range = [final_lons.min(), final_lons.max()]
lat_range = [final_lats.min(), final_lats.max()]

lon_bins = np.linspace(lon_range[0], lon_range[1], 150)
lat_bins = np.linspace(lat_range[0], lat_range[1], 150)

H, xedges, yedges = np.histogram2d(final_lons, final_lats, bins=[lon_bins, lat_bins])
H_smooth = gaussian_filter(H.T, sigma=3)

# 绘制
im = ax.imshow(
    H_smooth,
    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    origin='lower',
    cmap='hot',
    aspect='auto',
    alpha=0.8
)

# 标记排放口
ax.scatter(
    sources['longitude'],
    sources['latitude'],
    c='cyan',
    s=200,
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
        xytext=(5, 5),
        textcoords='offset points',
        fontsize=10,
        color='white',
        weight='bold',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='blue', alpha=0.8)
    )

plt.colorbar(im, ax=ax, label='Particle Density')
ax.set_xlabel('Longitude (°E)', fontsize=12)
ax.set_ylabel('Latitude (°N)', fontsize=12)
ax.set_title(
    f'Particle Density Distribution\n{START_DATE.strftime("%Y-%m-%d")} to {END_DATE.strftime("%Y-%m-%d")}',
    fontsize=14,
    weight='bold'
)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

density_file = os.path.join(OUTPUT_DIR, 'particle_density.png')
plt.savefig(density_file, dpi=300, bbox_inches='tight')
print(f"[OK] 密度图已保存: {density_file}")

# 轨迹图（采样）
plt.figure(figsize=(14, 12))
n_plot = min(50, n_particles)
sample_indices = np.random.choice(n_particles, n_plot, replace=False)

for idx in sample_indices:
    traj = trajectory[:, idx, :]
    valid = ~(np.isnan(traj[:, 0]) | np.isnan(traj[:, 1]))
    plt.plot(traj[valid, 0], traj[valid, 1], alpha=0.4, linewidth=0.8)

plt.scatter(
    sources['longitude'],
    sources['latitude'],
    c='red',
    s=250,
    marker='*',
    edgecolors='white',
    linewidths=2,
    label='Sources',
    zorder=10
)

plt.xlabel('Longitude (°E)', fontsize=12)
plt.ylabel('Latitude (°N)', fontsize=12)
plt.title(
    f'Particle Trajectories (Sample: {n_plot} particles)\n{START_DATE.strftime("%Y-%m-%d")} to {END_DATE.strftime("%Y-%m-%d")}',
    fontsize=14,
    weight='bold'
)
plt.legend(fontsize=10)
plt.grid(True, alpha=0.3)

traj_file = os.path.join(OUTPUT_DIR, 'particle_trajectories.png')
plt.savefig(traj_file, dpi=300, bbox_inches='tight')
print(f"[OK] 轨迹图已保存: {traj_file}")

# 关闭文件
ds.close()

print("\n" + "=" * 60)
print("完成！")
print("=" * 60)
print(f"结果保存在: {OUTPUT_DIR}")
print(f"- 轨迹数据: trajectory.pkl")
print(f"- 密度图: particle_density.png")
print(f"- 轨迹图: particle_trajectories.png")
print("=" * 60)
