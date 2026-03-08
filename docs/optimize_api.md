# `gbeampro.optimize` API Reference

## Table of Contents

- [Overview](#overview)
- [データクラス](#データクラス)
  - [`Operand`](#operand)
  - [`OptResult`](#optresult)
- [関数](#関数)
  - [`optimize_astigmatic`](#optimize_astigmatic)
  - [`waist_operands`](#waist_operands)
  - [`find_minimum_system`](#find_minimum_system)
  - [`build_xy_systems`](#build_xy_systems)
  - [`find_lens_system`](#find_lens_system)
- [使用例](#使用例)
  - [ウェストへの収束](#ウェストへの収束)
  - [最小レンズ枚数の探索](#最小レンズ枚数の探索)
  - [結果の可視化](#結果の可視化)

---

## Overview

Zemax OpticStudio の Merit Function / Optimization を参考にした，非点収差ガウシアンビーム系の最適化モジュール。

- **Merit function** をオペランドのリストとして記述
- **アルゴリズム** を `'de'` / `'lm'` / `'hammer'` から選択
- **非点収差系**（球面レンズ + シリンドリカルレンズ混在）に対応
- 離散焦点距離・レンズ間隔制約・焦点距離絶対値制限をサポート

---

## データクラス

### `Operand`

Merit function の 1 項を表す。

```python
@dataclass
class Operand:
    type:   str    # オペランド種別（下表参照）
    z_mm:   float  # 評価位置 (mm)
    target: float  # 目標値
    weight: float = 1.0  # 相対重み
```

**`type` 一覧**

| `type` | 評価量 | 単位 |
|--------|--------|------|
| `'wx'` | x 軸ビーム半径 | mm |
| `'wy'` | y 軸ビーム半径 | mm |
| `'cvx'` | x 軸波面曲率 1/R（`target=0` → ウェスト） | mm⁻¹ |
| `'cvy'` | y 軸波面曲率 1/R | mm⁻¹ |
| `'thx'` | x 軸半発散角 | mrad |
| `'thy'` | y 軸半発散角 | mrad |

**残差の正規化**

| type | スケール |
|------|---------|
| `wx`, `wy` | `target` |
| `cvx`, `cvy` | `1 / z_mm`（残差 ≈ z/R ≈ Δz/z） |
| `thx`, `thy` | `target` |

Merit value = Σ `weight × (val − target)² / scale²`

---

### `OptResult`

最適化結果。

```python
@dataclass
class OptResult:
    specs:          list[dict]   # レンズ仕様リスト（下記参照）
    merit:          float        # 最終 merit 値
    operand_values: list[float]  # 各オペランドの実測値
    success:        bool         # アルゴリズム収束フラグ
    elapsed_s:      float        # 実行時間 (s)
```

`specs` の各要素:

```python
{'type': 'spherical' | 'cyl_x' | 'cyl_y', 'z_mm': float, 'f_mm': float}
```

---

## 関数

### `optimize_astigmatic`

非点収差ガウシアンビーム系を最適化する主関数。

```python
optimize_astigmatic(
    beam: GaussBeam,
    lens_types: list[str],
    operands: list[Operand],
    *,
    z_bounds: tuple[float, float] | None = None,
    z_max_mm: float | None = None,
    f_bounds: tuple[float, float] | None = None,
    f_abs_bounds: tuple[float, float] | None = None,
    f_step_mm: float = 0.0,
    min_lens_sep_mm: float = 0.0,
    algorithm: str = 'de',
    seed: int = 42,
    maxiter: int = 1000,
    popsize: int = 15,
    tol: float = 1e-12,
) -> OptResult
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `beam` | 入力ビーム（x/y 共通） |
| `lens_types` | レンズ種別リスト。`'spherical'`・`'cyl_x'`・`'cyl_y'` |
| `operands` | Merit function オペランドのリスト |
| `z_bounds` | レンズ z 位置の探索範囲 (mm)。デフォルト: `(beam.z_mm, max operand z)` |
| `z_max_mm` | レンズ z 位置の上限 (mm)。`z_bounds` を更に絞る |
| `f_bounds` | 焦点距離の探索範囲 (mm)。デフォルト: `(1, 10000)` |
| `f_abs_bounds` | `|f|` の許容範囲 `(f_min, f_max)`。例: `(30, 1000)` → f ∈ [−1000,−30]∪[30,1000] |
| `f_step_mm` | 焦点距離の離散ステップ (mm)。`0` で連続 |
| `min_lens_sep_mm` | 隣接レンズ間の最小間隔 (mm) |
| `algorithm` | 最適化アルゴリズム（下表） |
| `seed` | 乱数シード |
| `maxiter` | 最大イテレーション数 |
| `popsize` | DE の個体数（DE / Hammer のみ） |
| `tol` | 収束判定値 |

**アルゴリズム**

| `algorithm` | 説明 |
|-------------|------|
| `'de'` | Differential Evolution — グローバル探索，ロバスト |
| `'lm'` | Levenberg-Marquardt / TRF — ローカル探索，高速 |
| `'hammer'` | DE グローバル → TRF 精密化（Zemax Hammer 相当） |

**レンズ種別と軸への作用**

| `type` | x 軸 | y 軸 |
|--------|------|------|
| `'spherical'` | ThinLens(f) | ThinLens(f) |
| `'cyl_x'` | ThinLens(f) | 透過（恒等） |
| `'cyl_y'` | 透過（恒等） | ThinLens(f) |

**実行時間**: 結果に `elapsed_s` として格納され，標準出力にも表示される。

---

### `waist_operands`

ウェスト収束を目標とする標準オペランドセットを生成する。

```python
waist_operands(
    z_mm: float,
    wx_mm: float,
    wy_mm: float,
    size_weight: float = 1.0,
    waist_tol_x_mm: float | None = None,
    waist_tol_y_mm: float | None = None,
    curvature_weight_x: float | None = None,
    curvature_weight_y: float | None = None,
) -> list[Operand]
```

生成されるオペランド: `wx`, `wy`, `cvx`（target=0）, `cvy`（target=0）

**曲率重み (`cvx`/`cvy`) の決定ロジック**

優先順位: `curvature_weight_*` 明示 > `waist_tol_*_mm` 換算 > `size_weight` と同値

`waist_tol_*_mm` を指定した場合:

```
cw = size_weight × (w_mm / (tol_mm / z_mm))²
```

これにより，ウェスト位置ずれ `tol_mm` が size residual と同程度の寄与になる。

> **注意**: `waist_tol_mm` が大きいと `cw` が極端に小さくなり，ウェスト位置の制約が事実上なくなる。`waist_tol_mm` が小さすぎると `cw` が巨大になり最適化が失敗しやすい。
> 不明な場合は `curvature_weight_x=1.0, curvature_weight_y=1.0` を明示的に指定することを推奨。

---

### `find_minimum_system`

Merit 閾値を満たす最小レンズ枚数を自動探索する。

```python
find_minimum_system(
    beam: GaussBeam,
    operands: list[Operand],
    lens_type: str | list[str] = 'spherical',
    max_lenses: int = 5,
    merit_threshold: float = 1e-3,
    # optimize_astigmatic と同じ制約引数...
    verbose: bool = True,
) -> tuple[OptResult, int]
```

`n = 1, 2, ...` と枚数を増やしながら `optimize_astigmatic` を呼び出し，
`merit < merit_threshold` を最初に満たした時点で停止する。

**戻り値**: `(OptResult, n_lenses)`

---

### `build_xy_systems`

レンズ仕様リストから x/y 軸の `OpticalSystem` を構築する。

```python
build_xy_systems(
    beam: GaussBeam,
    lens_specs: list[dict],
    z_target_mm: float,
) -> tuple[OpticalSystem, OpticalSystem]
```

`optimize_astigmatic` の結果を可視化・検証するために使う。

```python
sx, sy = build_xy_systems(beam, result.specs, z_target_mm + 100)
traj_x = sx.trace(beam, dz=1.0)
traj_y = sy.trace(beam, dz=1.0)
```

---

### `find_lens_system`

回転対称レンズ系の 1 軸最適化（シンプルな用途向け）。

```python
find_lens_system(
    beam: GaussBeam,
    z_target_mm: float,
    w_target_mm: float,
    n_lenses: int = 1,
    z_bounds: tuple[float, float] | None = None,
    f_bounds: tuple[float, float] = (1.0, 1000.0),
    tol: float = 1e-10,
    seed: int = 42,
) -> OpticalSystem
```

---

## 使用例

### ウェストへの収束

球面レンズ1枚＋シリンドリカルレンズ2枚で，円形ビームを楕円ウェストに収束させる例。

```python
from gbeampro import GaussBeam
from gbeampro.optimize import waist_operands, optimize_astigmatic, build_xy_systems

# 入力ビーム: λ=800 nm, ウェスト半径 w₀=1.0 mm
beam = GaussBeam.from_waist(wl_um=0.8, w0_mm=1.0)

# メリット関数: z=500 mm で wx=120 µm, wy=400 µm のウェストを目標
# curvature_weight はウェスト条件 (cvx=cvy=0) の重み
operands = waist_operands(
    z_mm=500, wx_mm=0.12, wy_mm=0.40,
    size_weight=1.0,
    curvature_weight_x=1.0,   # x 軸ウェスト条件の重み
    curvature_weight_y=1.0,   # y 軸ウェスト条件の重み
)

result = optimize_astigmatic(
    beam,
    lens_types=['spherical', 'cyl_x', 'cyl_y'],  # 球面 + シリンドリカル x + シリンドリカル y
    operands=operands,
    f_abs_bounds=(30, 1000),   # |f| を [30, 1000] mm に制限（市販カタログ相当）
    f_step_mm=5.0,             # 焦点距離を 5 mm 刻みに離散化
    z_max_mm=480.0,            # レンズ配置の上限位置（ターゲット手前 20 mm）
    min_lens_sep_mm=20.0,      # 隣接レンズ間の最小間隔
    algorithm='de',            # Differential Evolution によるグローバル探索
    maxiter=2000, popsize=20, seed=42,
)

print(result.specs)    # [{'type': ..., 'z_mm': ..., 'f_mm': ...}, ...]
print(f'merit={result.merit:.2e}  elapsed={result.elapsed_s:.1f}s')
```

### 最小レンズ枚数の探索

`n=1, 2, ...` と枚数を増やしながら最適化を繰り返し，merit が閾値を下回った時点で停止する。

```python
from gbeampro.optimize import find_minimum_system

result, n = find_minimum_system(
    beam, operands,
    lens_type='spherical',     # 同種レンズを n 枚使用
    max_lenses=5,              # 最大試行枚数
    merit_threshold=1e-3,      # この値を下回れば合格
    f_abs_bounds=(30, 1000),   # |f| を [30, 1000] mm に制限
    f_step_mm=5.0,             # 焦点距離を 5 mm 刻みに離散化
    z_max_mm=140.0,            # レンズ配置の上限位置
    min_lens_sep_mm=20.0,      # 隣接レンズ間の最小間隔
    algorithm='de',
)
print(f'最小枚数: {n} 枚,  merit={result.merit:.2e}')
```

### 結果の可視化

x/y 各軸の caustic を同一グラフに重ね描きする。

```python
import gbeampro.plot as gplot
import matplotlib.pyplot as plt

# result.specs から x/y 軸の光学系を構築
sx, sy = build_xy_systems(beam, result.specs, 600)
fig, ax = plt.subplots(figsize=(12, 4))
gplot.plot_system(sx, sx.trace(beam, dz=1.0), ax, label='x')
gplot.plot_system(sy, sy.trace(beam, dz=1.0), ax, label='y')
ax.axvline(500, color='k', ls=':')  # ターゲット位置を点線で表示
plt.tight_layout()
```
