# `gbeampro.analysis` / `gbeampro.plot` API Reference

## Table of Contents

- [Overview](#overview)
- [`gbeampro.analysis`](#gbeamproanalysis)
  - [`find_waists`](#find_waists)
  - [`rayleigh_range`](#rayleigh_range)
  - [`confocal_parameter`](#confocal_parameter)
  - [`beam_at`](#beam_at)
  - [`aperture_loss`](#aperture_loss)
- [`gbeampro.plot`](#gbeamproplot)
  - [`plot_caustic`](#plot_caustic)
  - [`plot_system`](#plot_system)
    - [素子シンボルの説明](#素子シンボルの説明)
- [使用例](#使用例)
  - [ウェスト検出と Rayleigh 長の計算](#ウェスト検出と-rayleigh-長の計算)
  - [アパーチャによるパワー損失](#アパーチャによるパワー損失)
  - [caustic の単純プロット](#caustic-の単純プロット)
  - [光学系の可視化（素子シンボル付き）](#光学系の可視化素子シンボル付き)
  - [複数ビームの重ね描き](#複数ビームの重ね描き)

---

## Overview

`gbeampro.analysis` はビーム軌跡の解析関数を提供し，`gbeampro.plot` は matplotlib を用いた可視化関数を提供する。

- `analysis` は `system.trace()` の戻り値（`list[GaussBeam]`）を入力として受け取る
- `plot` 関数は既存の `matplotlib.axes.Axes` オブジェクトに描画するため，柔軟なレイアウトに対応できる

---

## `gbeampro.analysis`

### `find_waists`

ビーム軌跡の中からウェスト位置を検出する。

```python
def find_waists(trajectory: list[GaussBeam]) -> list[GaussBeam]
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `trajectory` | `OpticalSystem.trace()` の戻り値 |

**戻り値**: ウェスト点の `GaussBeam` リスト

**検出ロジック**

波面曲率半径 R が負（収束）から正（発散）に転じる点をウェストとして検出する。ウェスト点では R が −∞ から +∞ へ遷移する（符号が反転する）。

> **注意**: 入力軌跡の `dz` が粗い場合，ウェスト位置の精度が低下する。精密な位置が必要な場合は `dz` を小さくすること。

---

### `rayleigh_range`

ビームの Rayleigh 長を計算する。

```python
def rayleigh_range(beam: GaussBeam) -> float
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `beam` | Rayleigh 長を計算する `GaussBeam`。ウェスト位置で呼び出すと w₀ に基づく値が得られる |

**戻り値**: Rayleigh 長 z_R (mm)

**計算式**

$$z_R = \frac{\pi \, n \, w^2}{\lambda}$$

w はビーム半径，λ は波長，n は屈折率。

---

### `confocal_parameter`

ビームの共焦点パラメータを計算する。

```python
def confocal_parameter(beam: GaussBeam) -> float
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `beam` | 計算対象の `GaussBeam` |

**戻り値**: 共焦点パラメータ b (mm)

**計算式**

$$b = 2 z_R = \frac{2 \pi \, n \, w^2}{\lambda}$$

共焦点パラメータは，ビームが焦点付近で回折限界以内に収まる軸方向の長さの指標となる。

---

### `beam_at`

軌跡上の任意の位置 z におけるビームを返す。

```python
def beam_at(trajectory: list[GaussBeam], z_mm: float) -> GaussBeam
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `trajectory` | `OpticalSystem.trace()` の戻り値 |
| `z_mm` | ビームを求める位置 (mm) |

**戻り値**: 位置 `z_mm` における `GaussBeam`

`z_mm` 以下で最後の軌跡点から `Propagation` で残り距離だけ自由伝搬させて求めるため，`dz` の刻みに依存しない。

- `z_mm` が軌跡の終端より先の場合は，自由空間として外挿する
- `z_mm` が軌跡の開始位置より前の場合は `ValueError`
- `z_mm` に薄肉素子（レンズ等）がある場合は素子通過**後**のビームを返す

---

### `aperture_loss`

光軸中心に置いた半径 r の円形アパーチャで遮られるパワーの割合を計算する。

```python
def aperture_loss(
    beam: GaussBeam,
    r_mm: float,
    beam_y: GaussBeam | None = None,
) -> float
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `beam` | アパーチャ位置でのビーム（x 方向）。`beam_at` で取得する |
| `r_mm` | アパーチャ半径 (mm) |
| `beam_y` | アパーチャ位置での y 方向ビーム。`None` の場合は円形ビーム（`w_y = w_x`）として扱う |

**戻り値**: 損失割合 L（0–1）。% 表示は `L * 100`

**計算式**

円形ビーム（w_x = w_y = w）:

$$L = \exp\left(-\frac{2r^2}{w^2}\right)$$

楕円ビーム（w_x ≠ w_y）: 強度 $I \propto \exp(-2x^2/w_x^2 - 2y^2/w_y^2)$ を極座標で表し，動径方向を解析積分する。

$$1 - L = \frac{1}{2\pi\, w_x w_y} \int_0^{2\pi} \frac{1 - \exp\left(-2r^2 a(\phi)\right)}{a(\phi)}\, d\phi, \qquad a(\phi) = \frac{\cos^2\phi}{w_x^2} + \frac{\sin^2\phi}{w_y^2}$$

角度方向の積分は周期関数の台形則（4096 点）で評価する。

> **注意**: アパーチャによる回折（下流のビーム形状の変化）は考慮しない。損失割合のみを返す。

---

## `gbeampro.plot`

### `plot_caustic`

ビーム軌跡から caustic 曲線（z に対するビーム半径 w の変化）をプロットする。

```python
def plot_caustic(
    trajectory: list[GaussBeam],
    ax,
    **kwargs,
)
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `trajectory` | `OpticalSystem.trace()` の戻り値 |
| `ax` | 描画先の `matplotlib.axes.Axes` |
| `**kwargs` | `ax.plot()` に渡す追加キーワード引数（`color`, `lw`, `label` など） |

軌跡の `(z_mm, w_mm)` を折れ線でプロットする。`+w` と `−w` の両側は描画されない（片側のみ）。

---

### `plot_system`

caustic と光学素子のシンボルを重ねてプロットする。複数回呼び出すと色が自動的に変わり，複数ビームの比較が容易になる。

```python
def plot_system(
    system,
    trajectory: list[GaussBeam],
    ax,
    beam_kw: dict | None = None,
    label: str = "",
)
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `system` | `OpticalSystem` インスタンス |
| `trajectory` | `system.trace()` の戻り値 |
| `ax` | 描画先の `matplotlib.axes.Axes` |
| `beam_kw` | caustic の描画オプション（`dict`）。`None` の場合はデフォルトスタイル |
| `label` | 凡例に表示するラベル文字列 |

`plot_caustic` を内部で呼び出してビーム半径を描画した後，各光学素子のシンボルを重ねて描画する。

#### 素子シンボルの説明

| 素子クラス | シンボル |
|-----------|---------|
| `ThinLens` | 矢印付き縦線（集光レンズ: ↕，発散レンズ: ↔ 相当） |
| `Interface` | 破線の縦線 |
| `InterfaceCurved` | 破線の縦線（`Interface` と同様） |
| `CurvedMirrorTan` / `CurvedMirrorSag` | 実線の縦線 |
| `Propagation` | シンボルなし（ビーム曲線のみ） |

> **複数ビームの重ね描き**: 同じ `ax` に対して `plot_system` を複数回呼び出すと，各ビームの caustic 色が matplotlib のデフォルトカラーサイクルに従って自動的に変わる。

---

## 使用例

### ウェスト検出と Rayleigh 長の計算

```python
from gbeampro import GaussBeam
from gbeampro.elements import Propagation, ThinLens
from gbeampro.system import OpticalSystem
from gbeampro.analysis import find_waists, rayleigh_range, confocal_parameter

beam = GaussBeam.from_waist(wl_um=1.064, w0_mm=1.0)

system = (OpticalSystem()
    .add(Propagation(200))
    .add(ThinLens(f_mm=150))
    .add(Propagation(300))
)

trajectory = system.trace(beam, dz=0.5)

# ウェストを検出
waists = find_waists(trajectory)
for w in waists:
    zR = rayleigh_range(w)
    b  = confocal_parameter(w)
    print(f'ウェスト: z={w.z_mm:.2f} mm,  w₀={w.w_mm:.4f} mm')
    print(f'  Rayleigh 長   z_R = {zR:.2f} mm')
    print(f'  共焦点パラメータ b = {b:.2f} mm')
```

### アパーチャによるパワー損失

z = 250 mm に半径 0.3 mm のアパーチャがある場合の損失を計算する。

```python
from gbeampro import GaussBeam
from gbeampro.optimize import build_xy_systems
from gbeampro.analysis import beam_at, aperture_loss

beam = GaussBeam.from_waist(wl_um=1.064, w0_mm=2.0)
specs = [
    {'type': 'spherical', 'z_mm': 0.77,   'f_mm': 140.0},
    {'type': 'spherical', 'z_mm': 84.55,  'f_mm': -245.0},
    {'type': 'spherical', 'z_mm': 136.16, 'f_mm': -30.0},
]

Z_APERTURE, R_APERTURE = 250.0, 0.3
sx, sy = build_xy_systems(beam, specs, Z_APERTURE)
bx = beam_at(sx.trace(beam, dz=0.5), Z_APERTURE)
by = beam_at(sy.trace(beam, dz=0.5), Z_APERTURE)

loss = aperture_loss(bx, R_APERTURE, beam_y=by)
print(f'w = {bx.w_mm*1e3:.1f} µm,  loss = {loss*100:.3f} %')
```

### caustic の単純プロット

```python
import matplotlib.pyplot as plt
from gbeampro import GaussBeam
from gbeampro.elements import Propagation, ThinLens
from gbeampro.system import OpticalSystem
import gbeampro.plot as gplot

beam = GaussBeam.from_waist(wl_um=0.8, w0_mm=0.5)

system = (OpticalSystem()
    .add(Propagation(100))
    .add(ThinLens(f_mm=80))
    .add(Propagation(150))
)

trajectory = system.trace(beam, dz=1.0)

fig, ax = plt.subplots(figsize=(10, 4))
gplot.plot_caustic(trajectory, ax, color='steelblue', lw=1.5, label='ビーム半径')
ax.set_xlabel('z (mm)')
ax.set_ylabel('w (mm)')
ax.legend()
plt.tight_layout()
plt.show()
```

### 光学系の可視化（素子シンボル付き）

```python
import matplotlib.pyplot as plt
from gbeampro import GaussBeam
from gbeampro.elements import Propagation, ThinLens, Interface
from gbeampro.system import OpticalSystem
import gbeampro.plot as gplot

beam = GaussBeam.from_waist(wl_um=1.064, w0_mm=1.0)

system = (OpticalSystem()
    .add(Propagation(150))
    .add(Interface(n1=1.0, n2=1.5))   # ガラス入射
    .add(Propagation(20))
    .add(Interface(n1=1.5, n2=1.0))   # ガラス出射
    .add(Propagation(100))
    .add(ThinLens(f_mm=200))
    .add(Propagation(300))
)

trajectory = system.trace(beam, dz=1.0)

fig, ax = plt.subplots(figsize=(12, 4))
# caustic + 素子シンボルを同時に描画
gplot.plot_system(system, trajectory, ax, label='beam')
ax.set_xlabel('z (mm)')
ax.set_ylabel('w (mm)')
ax.legend()
plt.tight_layout()
plt.show()
```

### 複数ビームの重ね描き

x/y 軸で光学系が異なる非点収差系の可視化。

```python
import matplotlib.pyplot as plt
from gbeampro import GaussBeam
from gbeampro.elements import Propagation, ThinLens
from gbeampro.system import OpticalSystem
import gbeampro.plot as gplot

beam = GaussBeam.from_waist(wl_um=0.8, w0_mm=1.0)

# x 軸系
sx = (OpticalSystem()
    .add(Propagation(100))
    .add(ThinLens(f_mm=120))
    .add(Propagation(200))
)

# y 軸系（焦点距離が異なる）
sy = (OpticalSystem()
    .add(Propagation(100))
    .add(ThinLens(f_mm=80))
    .add(Propagation(200))
)

traj_x = sx.trace(beam, dz=1.0)
traj_y = sy.trace(beam, dz=1.0)

fig, ax = plt.subplots(figsize=(12, 4))
# 同じ ax に重ね描き → 色が自動的に変わる
gplot.plot_system(sx, traj_x, ax, label='x')
gplot.plot_system(sy, traj_y, ax, label='y')
ax.set_xlabel('z (mm)')
ax.set_ylabel('w (mm)')
ax.legend()
plt.tight_layout()
plt.show()
```
