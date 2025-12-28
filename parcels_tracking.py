"""
使用 Parcels 进行拉格朗日粒子追踪
优化用于处理大型 Delft3D FM 输出文件
"""

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from parcels import FieldSet, ParticleSet, JITParticle, AdvectionRK4, plotTrajectoriesFile
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import warnings
warnings.filterwarnings('ignore')

# ==================== 配置参数 ====================
# 文件路径
MAP_FILE = r'E:\HK_model\202206_runtime\202206_June_run.dsproj_data\HK-DFM\output\HK-DFM_map.nc'
SOURCE_FILE = r'E:\HK_model\202206_runtime\point.xlsx'
OUTPUT_DIR = r'E:\HK_model\202206_runtime\particle_tracking\parcels_results'

# 时间设置（基准时间: 2022-01-01）
START_DATE = datetime(2022, 6, 1)   # 数据实际开始时间
END_DATE = datetime(2022, 7, 11)    # 用户要求结束时间
BASE_DATE = datetime(2022, 1, 1)

# 粒子追踪参数
N_PARTICLES_PER_SOURCE = 200        # 每个排放口的粒子数（减少以提高速度）
OUTPUT_DT = timedelta(hours=6)      # 输出时间间隔（6小时）
RUNTIME = (END_DATE - START_DATE).total_seconds()  # 运行时长（秒）

print("=" * 60)
print("Parcels 粒子追踪系统")
print("=" * 60)
print(f"开始时间: {START_DATE}")
print(f"结束时间: {END_DATE}")
print(f"运行时长: {RUNTIME/86400:.1f} 天")
print(f"每个排放口粒子数: {N_PARTICLES_PER_SOURCE}")
print("=" * 60)

# ==================== 1. 加载流场数据 ====================
print("\n[1/5] 正在加载流场数据...")
print("提示: 由于文件较大(307GB)，此步骤可能需要几分钟...")

# Parcels 支持 UGRID 格式，但需要特殊配置
# 对于 Delft3D FM，我们需要手动指定变量和维度
filenames = {'U': MAP_FILE, 'V': MAP_FILE}
variables = {'U': 'mesh2d_ucx', 'V': 'mesh2d_ucy'}
dimensions = {
    'lon': 'mesh2d_face_x',
    'lat': 'mesh2d_face_y',
    'time': 'time'
}

# 创建 FieldSet（使用深度平均流速）
# 注意：Delft3D FM 的 3D 数据在 Parcels 中处理较复杂，我们先用表层或深度平均
try:
    fieldset = FieldSet.from_netcdf(
        filenames,
        variables,
        dimensions,
        mesh='flat',  # 使用平面网格
        allow_time_extrapolation=False,
        deferred_load=True  # 延迟加载，节省内存
    )
    print("✓ 流场数据加载成功")
except Exception as e:
    print(f"✗ 流场加载失败: {e}")
    print("\n尝试使用简化方法...")
    # 备用方案：使用自定义读取
    import netCDF4 as nc
    ds = nc.Dataset(MAP_FILE)
    print(f"  网格单元数: {len(ds.variables['mesh2d_face_x'][:])}")
    print(f"  时间步数: {len(ds.variables['time'][:])}")
    print(f"  层数: {ds.variables['mesh2d_ucx'].shape}")
    ds.close()
    raise

# ==================== 2. 加载排放口位置 ====================
print("\n[2/5] 正在加载排放口位置...")
sources = pd.read_excel(SOURCE_FILE)
print(f"✓ 排放口数量: {len(sources)}")
print(sources)

# ==================== 3. 创建粒子集 ====================
print("\n[3/5] 正在创建粒子集...")

# 计算开始时间相对于基准时间的秒数
start_seconds = (START_DATE - BASE_DATE).total_seconds()

# 为每个排放口创建粒子
lons = []
lats = []
times = []

for idx, source in sources.iterrows():
    # 在每个排放口周围随机分布粒子（约±1km）
    source_lons = source['longitude'] + np.random.normal(0, 0.01, N_PARTICLES_PER_SOURCE)
    source_lats = source['latitude'] + np.random.normal(0, 0.01, N_PARTICLES_PER_SOURCE)

    lons.extend(source_lons)
    lats.extend(source_lats)
    times.extend([start_seconds] * N_PARTICLES_PER_SOURCE)

# 创建 ParticleSet
pset = ParticleSet(
    fieldset=fieldset,
    pclass=JITParticle,
    lon=lons,
    lat=lats,
    time=times
)

total_particles = len(lons)
print(f"✓ 总粒子数: {total_particles}")

# ==================== 4. 执行粒子追踪 ====================
print(f"\n[4/5] 开始粒子追踪...")
print(f"预计运行时间: 根据系统性能，可能需要10-30分钟")
print("进度将实时显示...\n")

# 输出文件
import os
os.makedirs(OUTPUT_DIR, exist_ok=True)
output_file = os.path.join(OUTPUT_DIR, 'particle_trajectories.zarr')

# 执行追踪（使用4阶龙格库塔方法）
pset.execute(
    AdvectionRK4,
    runtime=RUNTIME,
    dt=timedelta(minutes=30),  # 时间步长30分钟
    output_file=pset.ParticleFile(name=output_file, outputdt=OUTPUT_DT),
    verbose_progress=True
)

print("\n✓ 粒子追踪完成！")

# ==================== 5. 生成密度估计图 ====================
print("\n[5/5] 正在生成密度估计图...")

# 读取轨迹数据
try:
    from parcels import ParticleFile
    pfile = ParticleFile(output_file)

    # 获取最后时刻的粒子位置用于密度图
    final_lons = pfile.lon[:, -1]
    final_lats = pfile.lat[:, -1]

    # 移除 NaN 值
    valid_mask = ~(np.isnan(final_lons) | np.isnan(final_lats))
    final_lons = final_lons[valid_mask]
    final_lats = final_lats[valid_mask]

    # 创建密度图
    fig, ax = plt.subplots(figsize=(12, 10))

    # 使用 2D 直方图创建密度图
    lon_bins = np.linspace(final_lons.min(), final_lons.max(), 100)
    lat_bins = np.linspace(final_lats.min(), final_lats.max(), 100)

    H, xedges, yedges = np.histogram2d(final_lons, final_lats, bins=[lon_bins, lat_bins])

    # 高斯平滑
    from scipy.ndimage import gaussian_filter
    H_smooth = gaussian_filter(H.T, sigma=2)

    # 绘制密度图
    im = ax.imshow(
        H_smooth,
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        origin='lower',
        cmap='hot',
        aspect='auto',
        alpha=0.8
    )

    # 添加排放口位置
    ax.scatter(
        sources['longitude'],
        sources['latitude'],
        c='blue',
        s=100,
        marker='*',
        edgecolors='white',
        linewidths=2,
        label='Emission Sources',
        zorder=10
    )

    # 标注排放口名称
    for idx, source in sources.iterrows():
        ax.annotate(
            source['Point'],
            (source['longitude'], source['latitude']),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=9,
            color='white',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='blue', alpha=0.7)
        )

    plt.colorbar(im, ax=ax, label='Particle Density')
    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Latitude (°N)')
    ax.set_title(f'Particle Density Distribution\n{START_DATE.date()} to {END_DATE.date()}')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 保存图片
    density_file = os.path.join(OUTPUT_DIR, 'particle_density.png')
    plt.savefig(density_file, dpi=300, bbox_inches='tight')
    print(f"✓ 密度图已保存: {density_file}")

    # 绘制轨迹图（采样部分粒子）
    plt.figure(figsize=(12, 10))
    n_plot = min(100, total_particles)  # 最多绘制100条轨迹
    sample_indices = np.random.choice(total_particles, n_plot, replace=False)

    for idx in sample_indices:
        traj_lon = pfile.lon[idx, :]
        traj_lat = pfile.lat[idx, :]
        valid = ~(np.isnan(traj_lon) | np.isnan(traj_lat))
        plt.plot(traj_lon[valid], traj_lat[valid], alpha=0.3, linewidth=0.5)

    plt.scatter(
        sources['longitude'],
        sources['latitude'],
        c='red',
        s=200,
        marker='*',
        edgecolors='white',
        linewidths=2,
        label='Sources',
        zorder=10
    )

    plt.xlabel('Longitude (°E)')
    plt.ylabel('Latitude (°N)')
    plt.title(f'Particle Trajectories (Sample: {n_plot} particles)\n{START_DATE.date()} to {END_DATE.date()}')
    plt.legend()
    plt.grid(True, alpha=0.3)

    traj_file = os.path.join(OUTPUT_DIR, 'particle_trajectories.png')
    plt.savefig(traj_file, dpi=300, bbox_inches='tight')
    print(f"✓ 轨迹图已保存: {traj_file}")

except Exception as e:
    print(f"✗ 可视化失败: {e}")
    print("但轨迹数据已成功保存！")

print("\n" + "=" * 60)
print("完成！")
print("=" * 60)
print(f"结果保存在: {OUTPUT_DIR}")
print(f"- 轨迹数据: {output_file}")
print(f"- 密度图: particle_density.png")
print(f"- 轨迹图: particle_trajectories.png")
print("=" * 60)
