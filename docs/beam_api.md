# `gbeampro.beam` API Reference

## Table of Contents

- [Overview](#overview)
- [データクラス](#データクラス)
  - [`GaussBeam`](#gaussbeam)
    - [フィールド](#フィールド)
    - [プロパティ](#プロパティ)
      - [`.q`](#q)
      - [`.theta`](#theta)
    - [コンストラクタ](#コンストラクタ)
      - [`from_waist`](#from_waist)
      - [`from_q`](#from_q)
- [使用例](#使用例)
  - [ウェストからビームを生成する](#ウェストからビームを生成する)
  - [q パラメータからビームを生成する](#q-パラメータからビームを生成する)
  - [ビームパラメータを確認する](#ビームパラメータを確認する)

---

## Overview

ガウシアンビームを表す不変データクラス `GaussBeam` を提供する。

- **q パラメータ**（複素ビームパラメータ）によりビームを完全に記述
- **frozen dataclass** であるため，作成後にフィールドを変更できない
- ABCD 行列変換と組み合わせて光学系のビーム追跡に使用する

---

## データクラス

### `GaussBeam`

ガウシアンビームの状態を表す frozen dataclass。

```python
@dataclass(frozen=True)
class GaussBeam:
    wl_um: float         # 波長 (µm)
    n: float = 1.0       # 媒質の屈折率
    z_mm: float = 0.0    # 光軸上の現在位置 (mm)
    R_mm: float = inf    # 波面曲率半径 (mm)。正: 発散, 負: 収束, inf: 平面波（ウェスト）
    w_mm: float = 1.0    # ビーム半径 (1/e² 強度半径) (mm)
```

#### フィールド

| フィールド | 型 | 説明 | デフォルト |
|-----------|-----|------|-----------|
| `wl_um` | `float` | 波長 (µm) | — |
| `n` | `float` | 媒質の屈折率 | `1.0` |
| `z_mm` | `float` | 光軸上の現在位置 (mm) | `0.0` |
| `R_mm` | `float` | 波面曲率半径 (mm)。`inf` はウェスト位置（平面波面） | `inf` |
| `w_mm` | `float` | ビーム半径 (1/e² 強度半径) (mm) | `1.0` |

> **注意**: `R_mm` の符号は，ビームが進行方向（+z）に向かって発散する場合に正，収束する場合に負。ウェスト位置では `R_mm = inf`。

---

#### プロパティ

##### `.q`

複素ビームパラメータ q を返す。

```python
@property
def q(self) -> complex
```

**計算式**

$$\frac{1}{q} = \frac{1}{R} - i \frac{\lambda}{\pi n w^2}$$

コードでは以下のように計算される：

```python
q = 1.0 / (1j * wl_um * 1e-3 / (π * n * w_mm²) + 1.0 / R_mm)
```

q パラメータは ABCD 行列変換に直接使用される：

$$q' = \frac{A q + B}{C q + D}$$

---

##### `.theta`

ビームの半発散角（遠視野発散角）を返す。

```python
@property
def theta(self) -> float
```

**計算式**（単位: rad）

$$\theta = \arctan\!\left(\frac{\lambda}{\pi n w}\right)$$

ウェスト半径 w₀ が小さいほど，また波長が長いほど発散角が大きくなる。

---

#### コンストラクタ

##### `from_waist`

ウェスト半径を指定してビームを生成する。ウェスト位置では波面曲率半径 R = ∞ となる。

```python
@classmethod
def from_waist(
    cls,
    wl_um: float,
    w0_mm: float,
    z_mm: float = 0.0,
    n: float = 1.0,
) -> GaussBeam
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `wl_um` | 波長 (µm) |
| `w0_mm` | ウェスト半径 (mm) |
| `z_mm` | ウェスト位置の z 座標 (mm) |
| `n` | 媒質の屈折率 |

内部では `R_mm=inf, w_mm=w0_mm` として `GaussBeam` を生成する。

---

##### `from_q`

複素 q パラメータからビームを生成する。測定または解析的に q が既知の場合に使用する。

```python
@classmethod
def from_q(
    cls,
    wl_um: float,
    n: float,
    q: complex,
    z_mm: float = 0.0,
) -> GaussBeam
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `wl_um` | 波長 (µm) |
| `n` | 媒質の屈折率 |
| `q` | 複素ビームパラメータ |
| `z_mm` | 現在の z 座標 (mm) |

**変換式**

$$\frac{1}{q} = \frac{1}{R} - i \frac{\lambda}{\pi n w^2}$$

- `Re(1/q) = 0` の場合は `R_mm = inf`（ウェスト）
- `w_mm` は虚部から逆算：$w = \sqrt{\lambda / (\pi n \cdot |\mathrm{Im}(1/q)|)}$

---

## 使用例

### ウェストからビームを生成する

```python
from gbeampro import GaussBeam

# λ=1064 nm, ウェスト半径 w₀=0.5 mm, 位置 z=0 mm
beam = GaussBeam.from_waist(wl_um=1.064, w0_mm=0.5)

print(beam.w_mm)   # 0.5
print(beam.R_mm)   # inf（ウェスト位置）
print(beam.theta)  # 半発散角 (rad)
```

### q パラメータからビームを生成する

```python
import numpy as np
from gbeampro import GaussBeam

# 曲率半径 R=200 mm, ビーム半径 w=1.0 mm のビームを q から生成
wl_um = 0.8
n = 1.0
w_mm = 1.0
R_mm = 200.0

# 1/q を組み立てて q に変換
q_inv = 1.0 / R_mm - 1j * wl_um * 1e-3 / (np.pi * n * w_mm**2)
q = 1.0 / q_inv

beam = GaussBeam.from_q(wl_um=wl_um, n=n, q=q, z_mm=50.0)
print(f'R={beam.R_mm:.1f} mm,  w={beam.w_mm:.4f} mm')
```

### ビームパラメータを確認する

```python
import numpy as np
from gbeampro import GaussBeam

beam = GaussBeam.from_waist(wl_um=0.8, w0_mm=0.2)

# Rayleigh 長の手計算（参考）
zR = np.pi * beam.n * beam.w_mm**2 / (beam.wl_um * 1e-3)
print(f'Rayleigh 長 z_R = {zR:.2f} mm')
print(f'半発散角 θ = {np.degrees(beam.theta):.4f} °')
print(f'q = {beam.q}')
```
