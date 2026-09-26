# 复合层板受热自由弯曲计算 API

基于 Python + FastAPI 的复合层板自由热弯曲求解服务。

## 计算模型

各层自下而上排列，设距底面高度为 `z`，第 `i` 层在高度 `z` 处的应力为：

```
σ_i(z) = E_i · (ε₀ + κ·z − α_i·ΔT)
```

其中 `ε₀` 为底面应变、`κ` 为曲率、`α_i` 为热膨胀系数、`ΔT` 为统一温差。
自由弯曲要求全截面轴力 `N` 与对底面的弯矩 `M` 同时为零：

```
N = Σ∫σ_i dz = 0        M = Σ∫σ_i·z dz = 0
```

由此得到关于 `ε₀`、`κ` 的 2×2 线性方程组（`t` 为层厚，`z_c` 为层中心高度）：

```
[ A  B ] [ ε₀ ]   [ N_T ]      A = ΣE·t          N_T = ΣE·α·ΔT·t
[ B  D ] [ κ  ] = [ M_T ]      B = ΣE·t·z_c      M_T = ΣE·α·ΔT·t·z_c
                               D = ΣE·(z_top³−z_bot³)/3
```

厚度与模量均为正时系数矩阵行列式恒正，解存在且唯一。服务返回 `ε₀`、`κ`
以及各层上、下表面的应力。

### 安全温差区间

ΔT = 0 时各层应力为零，且每个层表面的应力都是 ΔT 的齐次线性函数
`σ = k·ΔT`（斜率 `k` 由 ΔT = 1 的一次求解得到）。结合每层的拉伸许用应力
`σ_t` 与压缩许用应力 `σ_c`（均为正的有限值），每个表面对 ΔT 的约束为

```
−σ_c ≤ k·ΔT ≤ σ_t
```

`k ≠ 0` 时该式给出一个下界候选和一个上界候选；`k = 0`（应力不随温差
变化）不约束区间。所有候选中最紧的下界与上界即构成以零温差为基准、
升温与降温均不超限的最大连续温差区间（任一表面斜率非零时两侧均有界；
全部斜率为零时两侧均无界，以 `null` 返回）。相同边界按层序、下表面
优先确定控制面。

## 快速启动

### Docker Compose（推荐，一键启动）

```bash
docker compose up --build
```

### 仅 Dockerfile

```bash
docker build -t thermal-bending-api .
docker run --rm -p 8000:8000 thermal-bending-api
```

### 本地开发（不用 Docker）

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 配置

无需任何配置文件或环境变量即可运行。可选环境变量：无（端口映射在
`docker-compose.yml` 中通过 `ports: "8000:8000"` 配置，如需修改宿主机
端口，改为例如 `"8080:8000"` 即可）。

## 接口地址

服务启动后监听 `http://localhost:8000`：

| 方法 | 路径 | 说明 |
| ---- | ---- | ---- |
| POST | `http://localhost:8000/solve` | 求解自由热弯曲 |
| POST | `http://localhost:8000/safe-temperature-range` | 求安全温差区间 |
| GET  | `http://localhost:8000/health` | 健康检查 |
| GET  | `http://localhost:8000/docs` | Swagger 交互式文档 |

### 请求格式（POST /solve）

- `delta_t`：统一温差（有限数值，单位如 K）
- `layers`：自下而上排列的层数组，**至少 2 层**；每层：
  - `thickness`：厚度，必须 > 0
  - `modulus`：弹性模量，必须 > 0
  - `cte`：热膨胀系数，必须为有限数值（可正可负可为零）

单位需自洽（例如厚度 m、模量 Pa、系数 1/K、温差 K，则应力为 Pa、曲率为 1/m）。

```bash
curl -s -X POST http://localhost:8000/solve \
  -H "Content-Type: application/json" \
  -d '{
    "delta_t": 100.0,
    "layers": [
      {"thickness": 0.001, "modulus": 2.10e11, "cte": 1.2e-5},
      {"thickness": 0.002, "modulus": 7.0e10,  "cte": 2.3e-5}
    ]
  }'
```

### 响应示例

```json
{
  "bottom_strain": 0.0010377880184331795,
  "curvature": 0.5474654377880186,
  "layers": [
    {"index": 0, "z_bottom": 0.0,   "z_top": 0.001, "stress_bottom": -34064516.129032336, "stress_top": 80903225.80645156},
    {"index": 1, "z_bottom": 0.001, "z_top": 0.003, "stress_bottom": -50032258.06451614, "stress_top": 26612903.225806456}
  ]
}
```

- `bottom_strain`：底面应变 ε₀
- `curvature`：曲率 κ
- `layers[i].stress_bottom / stress_top`：第 i 层下、上表面应力（index 0 为最底层）

### 请求格式（POST /safe-temperature-range）

- `layers`：自下而上排列的层数组，**至少 2 层**；每层：
  - `thickness`：厚度，必须 > 0
  - `modulus`：弹性模量，必须 > 0
  - `cte`：热膨胀系数，必须为有限数值（可正可负可为零）
  - `tensile_allowable`：拉伸许用应力，必须为有限值且严格 > 0
  - `compressive_allowable`：压缩许用应力（取正值），必须为有限值且严格 > 0

许用值非法（≤ 0、NaN、Infinity 或缺失）时返回 422。单位需与 `/solve`
自洽（许用应力与模量同单位，则返回温差与 `delta_t` 同单位）。

```bash
curl -s -X POST http://localhost:8000/safe-temperature-range \
  -H "Content-Type: application/json" \
  -d '{
    "layers": [
      {"thickness": 0.001, "modulus": 2.10e11, "cte": 1.2e-5,
       "tensile_allowable": 2.5e8, "compressive_allowable": 1.5e8},
      {"thickness": 0.002, "modulus": 7.0e10,  "cte": 2.3e-5,
       "tensile_allowable": 1.2e8, "compressive_allowable": 8.0e7}
    ]
  }'
```

### 响应示例

```json
{
  "delta_t_min": -185.40669856459323,
  "delta_t_max": 159.89684074790463,
  "lower_control": {"layer_index": 0, "surface": "top", "stress_type": "compression"},
  "upper_control": {"layer_index": 1, "surface": "bottom", "stress_type": "compression"}
}
```

- `delta_t_min` / `delta_t_max`：以零温差为基准的最大连续安全温差区间
  下、上界；该侧无界时为 `null`（仅当所有层表面应力均不随温差变化，
  此时两侧同时为 `null`）
- `lower_control` / `upper_control`：有限边界的控制信息——控制层
  `layer_index`（0 为最底层）、控制表面 `surface`（`bottom`/`top`）与
  拉压类型 `stress_type`（`tension`/`compression`）；对应边界无界时为
  `null`。相同边界按层序、下表面优先选择控制面

### 错误处理

非法参数返回 HTTP 422，并指明出错字段与原因，例如：

```json
{"detail": [{"field": "layers.0.thickness", "message": "Input should be greater than 0"}]}
```

会被拒绝的情形包括：层数少于 2、厚度或模量 ≤ 0、拉伸/压缩许用应力
≤ 0、任何参数为 NaN/Infinity、缺字段、类型错误，以及输入量级导致的
数值溢出。

## 项目结构

```
app/
  main.py      FastAPI 应用与路由、错误处理
  schemas.py   请求/响应模型（Pydantic 校验）
  solver.py    纯标准库求解器（平衡方程组求解与应力计算）
tests/
  test_solver.py  求解器自检（含双金属片 Timoshenko 解析解对比）
Dockerfile
docker-compose.yml
requirements.txt
```

本地运行求解器自检（无需安装依赖）：

```bash
python3 tests/test_solver.py
```
