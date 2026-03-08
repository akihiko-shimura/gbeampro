# `gbeampro.analysis` / `gbeampro.plot` API Reference

## Table of Contents

- [Overview](#overview)
- [`gbeampro.analysis`](#gbeamproanalysis)
  - [`find_waists`](#find_waists)
  - [`rayleigh_range`](#rayleigh_range)
  - [`confocal_parameter`](#confocal_parameter)
- [`gbeampro.plot`](#gbeamproplot)
  - [`plot_caustic`](#plot_caustic)
  - [`plot_system`](#plot_system)
    - [素子シンボルの説明](#素子シンボルの説明)
- [使用例](#使用例)
  - [ウェスト検出と Rayleigh 長の計算](#ウェスト検出と-rayleigh-長の計算)
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
