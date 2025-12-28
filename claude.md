# 简单粒子追踪系统 (Simple Particle Tracker)

## 📋 项目概述

本项目是一个专门为处理大型 Delft3D Flexible Mesh (DFM) 水动力模型输出而设计的拉格朗日粒子追踪系统。能够高效处理数百GB级别的NetCDF文件，模拟污染物或示踪剂在河口和海洋环境中的输运扩散过程。

**适用场景:**
- 河口污染物扩散模拟
- 海洋示踪剂追踪
- 搜救路径预测
- 油污扩散分析

**技术特点:**
- ✅ 支持超大型NetCDF文件（测试过307GB）
- ✅ 优化内存使用，采用检查点机制防止数据丢失
- ✅ 正确处理坐标转换（米 ↔ 经纬度）
- ✅ 支持3D分层流场（可选表层/底层/深度平均流速）
- ✅ 平流-扩散耦合模型
- ✅ 自动生成带地理底图的可视化结果

---

## 📂 文件结构

```
particle_tracking/
├── optimized_tracking.py              # 核心追踪脚本（主程序）
├── visualize_with_basemap.py          # 可视化脚本（生成图表）
├── PRD.shp (及相关文件)               # 珠江三角洲底图（shapefile）
├── optimized_results/                 # 结果输出文件夹
│   ├── trajectory.pkl                 # 完整轨迹数据
│   ├── particle_density_with_basemap.png      # 密度分布图
│   └── particle_trajectories_with_basemap.png # 轨迹图
└── claude.md                          # 本文档
```

---

## 🔬 工作原理

### 1. 粒子追踪数学模型

每个粒子的位置更新遵循平流-扩散方程：

```
dx/dt = u(x,y,t) + √(2K) * dW/dt
dy/dt = v(x,y,t) + √(2K) * dW/dt
```

其中：
- `u, v`: 流速分量（从Delft3D FM输出读取）
- `K`: 扩散系数（默认0.5 m²/s）
- `dW/dt`: 白噪声（模拟湍流扩散）

### 2. 关键技术细节

#### 坐标转换 ⭐（重要）

Delft3D FM输出的流速单位是 **m/s**，而粒子位置使用 **经纬度（°）**，必须进行转换：

```python
# 地球曲率修正
lat_rad = np.radians(粒子纬度)
meters_per_deg_lon = 111320 * cos(lat_rad)  # 经度方向，随纬度变化
meters_per_deg_lat = 111320                  # 纬度方向，恒定

# 流速转换为经纬度变化率
dlon/dt = u / meters_per_deg_lon
dlat/dt = v / meters_per_deg_lat

# 扩散项也需要转换
diffusion_lon = random_meters / meters_per_deg_lon
diffusion_lat = random_meters / meters_per_deg_lat
```

**常见错误:** 直接将米单位的位移加到经纬度上，导致粒子坐标变成巨大数值（如-14469）

#### 3D流场处理

Delft3D FM输出维度为 `[time, faces, layers]`，需要选择使用哪一层：

```python
# 选项1: 表层流速
u = mesh2d_ucx[time, face, 0]

# 选项2: 底层流速
u = mesh2d_ucx[time, face, -1]

# 选项3: 深度平均流速（推荐）
u_all_layers = mesh2d_ucx[time, face, :]
u = np.mean(u_all_layers)  # 忽略masked值
```

#### KD-Tree空间索引

使用KD树快速查找粒子所在网格单元：

```python
from scipy.spatial import cKDTree

# 建立索引（只需一次）
kdtree = cKDTree(网格中心坐标)

# 查询（每个粒子每个时间步）
distance, index = kdtree.query(粒子坐标, k=1)
if distance < 0.05°:  # 约5.5km
    流速 = mesh2d_ucx[time, index, :]
```

---

## 📥 输入数据说明

### 1. Delft3D FM 水动力输出文件 (MAP_FILE)

**文件格式:** NetCDF4 (`.nc`)
**文件名示例:** `HK-DFM_map.nc`
**典型大小:** 几GB到数百GB

**重要提示:** 由于NC文件通常非常大（本案例为307GB），**不建议上传到GitHub**。用户需要自行准备Delft3D FM模型运行后的输出文件。

#### 文件内容说明

该文件必须包含以下变量（UGRID格式）：

**1) 网格信息**
```python
Variables:
  mesh2d_face_x(nFaces)    # 网格单元中心X坐标（经度）
    - 单位: degrees_east
    - 示例: 113.0 - 115.0°E

  mesh2d_face_y(nFaces)    # 网格单元中心Y坐标（纬度）
    - 单位: degrees_north
    - 示例: 21.5 - 23.5°N

  mesh2d_node_x(nNodes)    # 网格节点X坐标
  mesh2d_node_y(nNodes)    # 网格节点Y坐标
  mesh2d_face_nodes(nFaces, nMaxNodes)  # 单元-节点连接关系
```

**2) 流速场（核心数据）**
```python
Variables:
  mesh2d_ucx(time, nFaces, nLayers)  # X方向流速分量
    - 单位: m/s
    - 维度: [时间, 网格单元, 垂向层数]
    - 示例形状: (2182, 101447, 10)
    - 数据类型: masked array (陆地为masked)

  mesh2d_ucy(time, nFaces, nLayers)  # Y方向流速分量
    - 单位: m/s
    - 维度: 同上
```

**3) 时间维度**
```python
Variables:
  time(time)  # 时间序列
    - 单位: seconds since BASE_DATE (如: 2022-01-01 00:00:00)
    - 示例: [0, 1800, 3600, ...] (每30分钟)
```

**4) 其他可选变量**
```python
Variables:
  mesh2d_waterdepth(time, nFaces)      # 水深
  mesh2d_s1(time, nFaces)              # 水位
  mesh2d_sa1(time, nFaces, nLayers)    # 盐度
  mesh2d_tem1(time, nFaces, nLayers)   # 温度
```

#### 如何查看NC文件内容

使用 `ncdump` 命令（NetCDF工具集）：
```bash
# 查看文件结构
ncdump -h your_model_map.nc

# 查看维度
ncdump -c your_model_map.nc | grep dimensions

# 查看变量列表
ncdump -v mesh2d_ucx your_model_map.nc | head
```

或使用Python：
```python
import netCDF4 as nc
ds = nc.Dataset('your_model_map.nc')
print(ds)  # 打印文件结构
print(ds.variables.keys())  # 列出所有变量
```

#### 数据要求检查清单

✅ 文件必须包含:
- [ ] `mesh2d_face_x` 和 `mesh2d_face_y`
- [ ] `mesh2d_ucx` 和 `mesh2d_ucy`（3D流速场）
- [ ] `time` 维度

✅ 数据格式:
- [ ] 坐标系统为经纬度（degrees）
- [ ] 流速单位为 m/s
- [ ] 时间单位为 seconds since BASE_DATE

❌ 不支持:
- 投影坐标系统（UTM等）需要先转换
- 2D流速场（需要至少1个垂向层）

#### 示例：本项目使用的文件

```
文件名: HK-DFM_map.nc
大小: 307 GB
网格单元数: 101,447
时间步数: 2,182 (每30分钟一个)
垂向层数: 10
时间范围: 2022-06-01 至 2022-07-16
覆盖区域: 珠江三角洲及邻近海域
  - 经度: 111.7°E - 115.6°E
  - 纬度: 21.5°N - 23.3°N
```

---

### 2. 排放口位置文件 (SOURCE_FILE)

**文件格式:** Excel (`.xlsx`) 或 CSV
**文件名示例:** `point.xlsx`

#### 文件结构

必须包含以下三列：

| 列名 | 数据类型 | 说明 | 示例 |
|------|----------|------|------|
| `Point` | 字符串 | 排放口名称 | "Humen", "Jiaomen" |
| `longitude` | 浮点数 | 经度（°E） | 113.4958 |
| `latitude` | 浮点数 | 纬度（°N） | 22.5755 |

#### 示例文件内容

```excel
Point      longitude   latitude
Humen      113.4958    22.5755
Jiaomen    113.4804    22.7935
Hongqili   113.3814    22.7625
Hengmen    113.4958    22.5755
Deep       114.0302    22.5091
```

#### 坐标要求
- 必须与NC文件使用相同的坐标系统（经纬度）
- 经纬度必须在NC文件覆盖范围内
- 精度建议保留4位小数（约11米精度）

---

### 3. 底图文件 (SHP_FILE, 可选)

**文件格式:** Shapefile (`.shp` 及相关文件)
**文件名示例:** `PRD.shp`

#### Shapefile组成

Shapefile是一组文件，必须包含：
```
PRD.shp     # 主文件（几何数据）
PRD.shx     # 索引文件
PRD.dbf     # 属性表
PRD.prj     # 投影信息
PRD.cpg     # 字符编码（可选）
PRD.sbn/sbx # 空间索引（可选）
```

#### 内容说明
- **几何类型:** 多边形（Polygon）或多部分多边形（MultiPolygon）
- **坐标系统:** WGS84（EPSG:4326）或与NC文件一致
- **覆盖范围:** 应包含粒子追踪区域

#### 如何准备底图
1. **下载:** GADM、Natural Earth、OpenStreetMap
2. **自定义:** 使用QGIS等GIS软件制作
3. **格式转换:** `ogr2ogr` 或 GeoPandas

示例（GeoPandas读取检查）：
```python
import geopandas as gpd
gdf = gpd.read_file('PRD.shp')
print(gdf.crs)  # 检查坐标系
print(gdf.total_bounds)  # 查看范围
gdf.plot()  # 预览
```

---

## 🚀 使用方法

### 环境要求

```bash
Python >= 3.8
numpy
pandas
scipy
netCDF4
matplotlib
geopandas
tqdm
```

安装依赖：
```bash
pip install numpy pandas scipy netCDF4 matplotlib geopandas tqdm
```

### 快速开始

#### 步骤0: 准备输入文件

确保你有：
1. ✅ Delft3D FM输出的NetCDF文件（`*_map.nc`）
2. ✅ 排放口位置Excel文件（`point.xlsx`）
3. ✅ （可选）区域底图Shapefile（`*.shp`）

#### 步骤1: 配置参数

编辑 `optimized_tracking.py` 文件中的配置部分（第21-37行）：

```python
# 文件路径
MAP_FILE = r'path/to/your/model_map.nc'      # Delft3D FM输出文件
SOURCE_FILE = r'path/to/emission_points.xlsx' # 排放口位置Excel文件
OUTPUT_DIR = r'path/to/output'                # 结果输出文件夹

# 时间范围（根据你的数据调整）
START_DATE = datetime(2022, 6, 1)
END_DATE = datetime(2022, 7, 11)
BASE_DATE = datetime(2022, 1, 1)  # NetCDF文件的时间基准

# 粒子参数
N_PARTICLES_PER_SOURCE = 100     # 每个排放口释放的粒子数
DT = 3600                        # 时间步长（秒），1小时
SAVE_INTERVAL = 48               # 保存间隔（步），每2天保存一次
DIFFUSION_COEF = 0.5             # 扩散系数（m²/s）
VELOCITY_LAYER = 'average'       # 'surface', 'bottom', 或 'average'
CHECKPOINT_INTERVAL = 200        # 检查点间隔（步）
```

#### 步骤2: 准备排放口文件

创建Excel文件（如 `point.xlsx`），包含以下列：

| Point    | longitude | latitude |
|----------|-----------|----------|
| Source1  | 113.4958  | 22.5755  |
| Source2  | 113.4804  | 22.7935  |
| ...      | ...       | ...      |

#### 步骤3: 运行追踪

```bash
python optimized_tracking.py
```

**预期输出:**
```
============================================================
优化版粒子追踪系统
============================================================

[1/4] 加载数据...
[OK] 网格单元数: 101447
[OK] 总时间步数: 2182
[OK] 使用时间步范围: 0 - 1920 (1921 步)

[2/4] 释放粒子...
[OK] 总粒子数: 500

[3/4] 开始粒子追踪...
追踪进度: 100%|██████████| 1921/1921 [01:38<00:00, 19.56it/s]

[OK] 追踪完成！实际保存了 41 个时间点
[OK] 结果已保存: trajectory.pkl

[4/4] 生成密度图...
[OK] 密度图已保存: particle_density.png
[OK] 轨迹图已保存: particle_trajectories.png
```

#### 步骤4: 生成带底图的可视化

```bash
python visualize_with_basemap.py
```

---

## ⚙️ 参数调优指南

### 内存和性能权衡

| 参数 | 增大 → | 减小 → |
|------|--------|--------|
| `N_PARTICLES_PER_SOURCE` | 更平滑的密度图，但内存↑ | 更快，内存↓ |
| `DT` | 更快，但精度↓ | 更精确，但时间↑ |
| `SAVE_INTERVAL` | 更少输出，内存↓ | 更详细时间序列 |
| `DIFFUSION_COEF` | 扩散更快 | 更集中跟随水流 |

### 大文件处理建议

**如果NetCDF文件 > 100GB:**
1. 减少粒子数至50-100/源
2. 增大时间步长至 3600-7200秒（1-2小时）
3. 增大保存间隔至 72-96步（3-4天）
4. 使用检查点机制（自动）

**如果内存溢出:**
```python
# 方法1: 减少同时读取的时间步数
# 已在代码中实现，每次只读取1个时间步

# 方法2: 使用分段运行
# 修改 START_DATE 和 END_DATE 分段模拟
```

---

## 📊 结果解读

### 输出文件说明

#### 1. `trajectory.pkl` (2-10 MB) ⭐ **重要**

这是最重要的输出文件，包含完整的粒子追踪结果。

**为什么重要？**
- ✅ **可重复分析**：无需重新运行模拟，直接加载数据进行各种分析
- ✅ **灵活可视化**：可以生成不同类型的图表（密度、轨迹、动画等）
- ✅ **数据共享**：文件小（~2MB），方便分享和存储
- ✅ **后续处理**：可以用于统计分析、路径计算、影响范围评估等

**文件格式：** Python pickle格式

**数据结构：**

```python
import pickle
with open('trajectory.pkl', 'rb') as f:
    data = pickle.load(f)

# 数据结构
data = {
    'trajectory': np.array,     # 形状: (n_times, n_particles, 2)
                                 # [时间点, 粒子ID, (经度, 纬度)]
                                 # 示例: (41, 500, 2) = 41个时间点，500个粒子

    'time_indices': np.array,   # 对应的NetCDF时间索引
                                 # 形状: (n_times,)
                                 # 用于回溯到原始NC文件

    'source_ids': np.array,     # 每个粒子的排放源ID
                                 # 形状: (n_particles,)
                                 # 示例: [0,0,0,...,1,1,1,...] (每个源100个粒子)

    'sources': pd.DataFrame,    # 排放源信息表
                                 # 列: Point, longitude, latitude

    'time': np.array,           # 完整时间序列（秒）
                                 # 从BASE_DATE开始计算

    'face_x': np.array,         # 网格中心x坐标（经度）
    'face_y': np.array          # 网格中心y坐标（纬度）
}
```

**使用示例：**

```python
import pickle
import numpy as np
import matplotlib.pyplot as plt

# 1. 加载数据
with open('trajectory.pkl', 'rb') as f:
    data = pickle.load(f)

trajectory = data['trajectory']  # (41, 500, 2)
sources = data['sources']
source_ids = data['source_ids']

# 2. 提取某个排放源的粒子
source_0_particles = np.where(source_ids == 0)[0]  # 第一个排放源的粒子索引
source_0_traj = trajectory[:, source_0_particles, :]  # 该源的轨迹

# 3. 计算最大扩散距离
from scipy.spatial.distance import cdist

for i, source in sources.iterrows():
    source_pos = np.array([[source['longitude'], source['latitude']]])
    final_positions = trajectory[-1, source_ids == i, :]  # 该源的最终位置

    # 计算距离（近似，未考虑地球曲率）
    distances = cdist(source_pos, final_positions, metric='euclidean')
    max_distance = distances.max()

    print(f"{source['Point']}: 最大扩散距离 {max_distance*111:.1f} km")

# 4. 计算粒子停留时间（在某个区域内）
region = {'lon': (113.0, 113.5), 'lat': (21.5, 22.0)}

for t in range(trajectory.shape[0]):
    pos_t = trajectory[t, :, :]
    in_region = (
        (pos_t[:, 0] >= region['lon'][0]) & (pos_t[:, 0] <= region['lon'][1]) &
        (pos_t[:, 1] >= region['lat'][0]) & (pos_t[:, 1] <= region['lat'][1])
    )
    print(f"时间步{t}: {in_region.sum()}个粒子在区域内")

# 5. 生成自定义可视化
plt.figure(figsize=(10, 8))
for i in range(0, 500, 10):  # 每10个粒子绘制一个
    traj = trajectory[:, i, :]
    plt.plot(traj[:, 0], traj[:, 1], alpha=0.5, linewidth=0.5)
plt.scatter(sources['longitude'], sources['latitude'],
            c='red', s=100, marker='*', zorder=10)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Custom Trajectory Visualization')
plt.savefig('custom_plot.png', dpi=200)
```

**建议上传到GitHub**：是的，文件小（~2MB），应该上传，方便用户：
- 直接查看示例结果
- 学习如何使用pkl文件
- 验证代码是否正确运行

#### 2. 密度分布图

显示模拟结束时刻粒子的空间分布：
- **热力图**: 红色区域表示高粒子密度
- **青色星号**: 排放源位置
- **灰色区域**: 陆地（底图）

**解读:**
- 高密度区 = 污染物累积区
- 从排放源到高密度区的路径 = 主要输运方向

#### 3. 轨迹图

显示采样粒子的完整运动轨迹：
- **彩色线条**: 每条线代表一个粒子的路径
- **往复运动**: 表示潮汐作用
- **定向输运**: 表示余流作用

---

## 🔧 常见问题

### Q1: 粒子向四周扩散，不跟随水流？

**原因:** 坐标转换错误

**解决:** 检查代码第160-174行，确保使用了正确的坐标转换：
```python
lat_rad = np.radians(pos[:, 1])
meters_per_deg_lon = 111320 * np.cos(lat_rad)
```

### Q2: 模拟运行到一半停止？

**原因:** 大文件读取超时或内存不足

**解决:**
1. 检查 `optimized_results/checkpoint_step_*.pkl` 文件
2. 减小 `N_PARTICLES_PER_SOURCE` 和增大 `DT`
3. 重新运行会自动从检查点恢复

### Q3: 密度图全是灰色/没有热力图？

**原因:** 粒子超出底图范围或坐标错误

**解决:**
1. 检查粒子坐标是否合理（经度113-115°, 纬度21-23°）
2. 打印 `final_lons` 和 `final_lats` 查看实际值
3. 调整 `visualize_with_basemap.py` 中的范围

### Q4: "distance too far" 导致粒子停止移动？

**原因:** 粒子离开了计算网格

**解决:**
```python
# 修改 optimized_tracking.py 第140行
if dist > 0.05:  # 减小这个阈值，或添加边界反射
    # 方案1: 忽略（当前）
    continue
    # 方案2: 边界反射（可选）
    # pos[i] = 上一个有效位置
```

---

## 🎓 技术背景

### 拉格朗日 vs 欧拉方法

| 方法 | 描述 | 优点 | 缺点 |
|------|------|------|------|
| **拉格朗日** | 追踪单个粒子 | 直观、适合示踪剂 | 计算量大 |
| **欧拉** | 固定网格浓度场 | 快速、适合浓度分布 | 数值扩散大 |

本项目使用拉格朗日方法。

### 相关工具对比

| 工具 | 优点 | 缺点 |
|------|------|------|
| **Parcels** | 功能强大，支持多种模型 | 对UGRID支持有限 |
| **Tracmass** | Fortran实现，速度快 | 配置复杂 |
| **本项目** | 简单、针对Delft3D FM优化 | 功能相对基础 |

---

## 📚 参考文献

1. Van Sebille, E., et al. (2018). "Lagrangian ocean analysis: Fundamentals and practices." *Ocean Modelling*, 121, 49-75.

2. Deltares (2023). "Delft3D-FM User Manual." Deltares, Delft, Netherlands.

3. Csanady, G.T. (1973). "Turbulent Diffusion in the Environment." D. Reidel Publishing Company.

---

## 📝 版本历史

### v1.0 (2024-12)
- ✅ 初始版本
- ✅ 支持Delft3D FM UGRID格式
- ✅ 坐标转换修复
- ✅ 3D流场处理
- ✅ 检查点机制
- ✅ 底图可视化

---

## 🤝 贡献

欢迎提交Issue和Pull Request！

**改进方向:**
- [ ] 支持更多水动力模型格式（ROMS, FVCOM等）
- [ ] GPU加速
- [ ] 粒子沉降和衰减
- [ ] 交互式可视化（Plotly/Bokeh）
- [ ] 并行计算（MPI）

---

## 📧 联系

如有问题或建议，请在GitHub Issue中讨论。

---

## ⚖️ 许可证

MIT License - 自由使用和修改

---

**最后更新:** 2024-12-28
**测试环境:** Windows 10, Python 3.13, Delft3D FM 1.2.95
