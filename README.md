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

由于 `ε₀`、`κ` 均为 `ΔT` 的线性函数且在 `ΔT = 0` 时为零，各层表面应力满足
`σ(ΔT) = m·ΔT`（斜率 `m` 即 `ΔT = 1` 时的应力）。给定每层拉伸/压缩许用
应力后，每个非零斜率表面解析给出正、负两个温差边界，零斜率表面不约束
区间；所有边界的交集即以零温差为基准、升温与降温均不超限的最大连续
温差区间（见 `POST /safe-temperature-range`）。

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

无需任何配置文件或环境变量即可运行。服务监听容器内 **8000** 端口，
`docker-compose.yml` 通过 `ports: "8000:8000"` 将其映射到宿主机 8000
端口；如需修改宿主机端口，改为例如 `"8080:8000"` 即可（容器内端口不变，
所有接口路径保持不变）。

## 接口地址

服务启动后监听 `http://localhost:8000`：

| 方法 | 路径 | 说明 |
| ---- | ---- | ---- |
| POST | `http://localhost:8000/solve` | 求解自由热弯曲 |
| POST | `http://localhost:8000/safe-temperature-range` | 求应力不超限的最大温差区间 |
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

求以零温差为基准、升温与降温均不超过许用应力的最大连续温差区间。

- `layers`：自下而上排列的层数组，**至少 2 层**；每层在 `thickness`、
  `modulus`、`cte`（要求同 `/solve`）之外还需：
  - `tensile_allowable`：拉伸许用应力，必须有限且 > 0
  - `compressive_allowable`：压缩许用应力（幅值），必须有限且 > 0，
    即应力不得小于 `-compressive_allowable`

```bash
curl -s -X POST http://localhost:8000/safe-temperature-range \
  -H "Content-Type: application/json" \
  -d '{
    "layers": [
      {"thickness": 0.001, "modulus": 2.10e11, "cte": 1.2e-5,
       "tensile_allowable": 2.5e8, "compressive_allowable": 2.5e8},
      {"thickness": 0.002, "modulus": 7.0e10,  "cte": 2.3e-5,
       "tensile_allowable": 1.5e8, "compressive_allowable": 1.5e8}
    ]
  }'
```

### 响应示例

```json
{
  "delta_t_min": -299.8065764023212,
  "delta_t_max": 299.8065764023212,
  "min_control": {"layer_index": 1, "surface": "bottom", "stress_type": "tension"},
  "max_control": {"layer_index": 1, "surface": "bottom", "stress_type": "compression"}
}
```

- `delta_t_min` / `delta_t_max`：区间下、上界；该方向无约束时为 `null`
  （例如各层热膨胀系数全部相同，任意温差下应力均为零）
- `min_control` / `max_control`：有限边界的控制位置（无界时为 `null`）：
  - `layer_index`：控制层序号（0 为最底层）
  - `surface`：控制表面，`bottom` 或 `top`
  - `stress_type`：边界处达到的限值类型，`tension`（拉）或 `compression`（压）
- 多个表面给出相同边界值时，按层序取序号最小者，同层下表面优先

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
  solver.py    纯标准库求解器（平衡方程组求解、应力计算与安全温差区间）
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
