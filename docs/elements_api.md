# `gbeampro.elements` API Reference

## Table of Contents

- [Overview](#overview)
- [抽象基底クラス](#抽象基底クラス)
  - [`Element`](#element)
    - [カスタム素子の作り方](#カスタム素子の作り方)
- [具体クラス](#具体クラス)
  - [`Propagation`](#propagation)
  - [`ThinLens`](#thinlens)
  - [`Interface`](#interface)
  - [`InterfaceCurved`](#interfacecurved)
  - [`CurvedMirrorTan`](#curvedmirrortan)
  - [`CurvedMirrorSag`](#curvedmirrorsag)
- [使用例](#使用例)
  - [基本的な素子の使用](#基本的な素子の使用)
  - [界面を通したビーム伝搬](#界面を通したビーム伝搬)
  - [カスタム素子の定義](#カスタム素子の定義)

---

## Overview

光学素子を ABCD 行列として表現するクラス群を提供する。

- すべての素子は抽象基底クラス `Element` を継承する
- `element.apply(beam)` により，ビームに ABCD 行列変換を適用できる
- `OpticalSystem.add()` と組み合わせて光学系を構築する

**ABCD 行列変換**

$$q' = \frac{A q + B}{C q + D}$$

---

## 抽象基底クラス

### `Element`

すべての光学素子の基底クラス。

```python
from abc import ABC, abstractmethod

class Element(ABC):
    @property
    @abstractmethod
    def matrix(self) -> np.ndarray: ...  # 2x2 ABCD 行列を返す

    def apply(self, beam: GaussBeam) -> GaussBeam:
        # q パラメータに ABCD 行列変換を適用して新しいビームを返す
        ...
```

**メソッド**

| メソッド | 説明 |
|---------|------|
| `matrix` | 2×2 ABCD 行列を `np.ndarray` で返す（抽象プロパティ） |
| `apply(beam)` | `beam` に ABCD 行列変換を適用し，新しい `GaussBeam` を返す |

`apply` は内部で `GaussBeam.from_q` を呼び出す。変換後の z 座標は `beam.z_mm + self.dz` に更新される。

---

#### カスタム素子の作り方

`Element` を継承し，`matrix` プロパティを実装するだけでカスタム素子を定義できる。

```python
import numpy as np
from gbeampro.elements import Element

class MyElement(Element):
    def __init__(self, param: float):
        self.param = param

    @property
    def matrix(self) -> np.ndarray:
        # ABCD 行列を返す
        A, B, C, D = 1.0, 0.0, -1.0 / self.param, 1.0
        return np.array([[A, B], [C, D]])
```

---

## 具体クラス

### `Propagation`

自由空間伝搬（距離 d_mm）。

```python
class Propagation(Element):
    def __init__(self, d_mm: float)
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `d_mm` | 伝搬距離 (mm) |

**ABCD 行列**

$$M = \begin{pmatrix} 1 & d \\ 0 & 1 \end{pmatrix}$$

ビーム半径と波面曲率が伝搬とともに変化する。

---

### `ThinLens`

薄肉レンズ（焦点距離 f_mm）。

```python
class ThinLens(Element):
    def __init__(self, f_mm: float)
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `f_mm` | 焦点距離 (mm)。正: 集光, 負: 発散 |

**ABCD 行列**

$$M = \begin{pmatrix} 1 & 0 \\ -1/f & 1 \end{pmatrix}$$

レンズはビーム半径を変えず，波面曲率のみを変化させる。

---

### `Interface`

平坦屈折界面（屈折率 n1 → n2）。

```python
class Interface(Element):
    def __init__(self, n1: float, n2: float)
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `n1` | 入射側の屈折率 |
| `n2` | 透過側の屈折率 |

**ABCD 行列**

$$M = \begin{pmatrix} 1 & 0 \\ 0 & n_1/n_2 \end{pmatrix}$$

変換後のビームの屈折率は n2 に更新される。ビーム半径は変化しないが，q パラメータが変わるため Rayleigh 長が変化する。

---

### `InterfaceCurved`

曲率を持つ屈折界面。

```python
class InterfaceCurved(Element):
    def __init__(self, n1: float, n2: float, r_mm: float)
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `n1` | 入射側の屈折率 |
| `n2` | 透過側の屈折率 |
| `r_mm` | 曲率半径 (mm)。正: 凸面（入射側から見て中心が遠い），負: 凹面 |

**ABCD 行列**

$$M = \begin{pmatrix} 1 & 0 \\ (n_2 - n_1)/(n_2 \cdot r) & n_1/n_2 \end{pmatrix}$$

平坦界面 `Interface` に集光・発散作用が加わった素子。

---

### `CurvedMirrorTan`

傾斜入射曲率ミラーのタンジェンシャル（tangential）面。

```python
class CurvedMirrorTan(Element):
    def __init__(self, r_mm: float, theta_deg: float)
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `r_mm` | 鏡面の曲率半径 (mm) |
| `theta_deg` | 入射角 (deg)。法線からの角度 |

**有効半径とABCD 行列**

タンジェンシャル面では入射角により有効半径が短縮される：

$$r_{\mathrm{eff}} = r \cdot \cos\theta$$

$$M = \begin{pmatrix} 1 & 0 \\ -2/r_{\mathrm{eff}} & 1 \end{pmatrix}$$

> **注意**: `CurvedMirrorTan` と `CurvedMirrorSag` を別々に使用することで，傾斜入射ミラーによる非点収差を正確にモデル化できる。

---

### `CurvedMirrorSag`

傾斜入射曲率ミラーのサジタル（sagittal）面。

```python
class CurvedMirrorSag(Element):
    def __init__(self, r_mm: float, theta_deg: float)
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `r_mm` | 鏡面の曲率半径 (mm) |
| `theta_deg` | 入射角 (deg)。法線からの角度 |

**有効半径とABCD 行列**

サジタル面では入射角により有効半径が伸長される：

$$r_{\mathrm{eff}} = r / \cos\theta$$

$$M = \begin{pmatrix} 1 & 0 \\ -2/r_{\mathrm{eff}} & 1 \end{pmatrix}$$

---

## 使用例

### 基本的な素子の使用

```python
from gbeampro import GaussBeam
from gbeampro.elements import Propagation, ThinLens

# λ=1064 nm, ウェスト半径 1.0 mm のビームを生成
beam = GaussBeam.from_waist(wl_um=1.064, w0_mm=1.0)

# 100 mm 伝搬後，f=150 mm のレンズを通過させる
p = Propagation(d_mm=100.0)
lens = ThinLens(f_mm=150.0)

beam1 = p.apply(beam)           # 伝搬後のビーム
beam2 = lens.apply(beam1)       # レンズ通過後のビーム

print(f'レンズ後: R={beam2.R_mm:.1f} mm,  w={beam2.w_mm:.4f} mm')
```

### 界面を通したビーム伝搬

```python
from gbeampro import GaussBeam
from gbeampro.elements import Propagation, Interface

beam = GaussBeam.from_waist(wl_um=1.064, w0_mm=0.5, n=1.0)

# 空気 (n=1) からガラス (n=1.5) への平坦界面
iface = Interface(n1=1.0, n2=1.5)
beam_in_glass = iface.apply(beam)

print(f'界面前: n={beam.n}, w={beam.w_mm:.4f} mm')
print(f'界面後: n={beam_in_glass.n}, w={beam_in_glass.w_mm:.4f} mm')

# ガラス内を 50 mm 伝搬
p = Propagation(d_mm=50.0)
beam_after = p.apply(beam_in_glass)
```

### カスタム素子の定義

```python
import numpy as np
from gbeampro import GaussBeam
from gbeampro.elements import Element

# 薄い波長板（位相のみ，ビームパラメータに影響しない恒等変換の例）
class IdentityElement(Element):
    @property
    def matrix(self) -> np.ndarray:
        return np.eye(2)  # 恒等 ABCD 行列

# 傾斜入射曲率ミラーで非点収差を生じさせる例
from gbeampro.elements import CurvedMirrorTan, CurvedMirrorSag

beam = GaussBeam.from_waist(wl_um=0.8, w0_mm=1.0)

# 曲率半径 500 mm, 入射角 10° のミラー（タンジェンシャル/サジタル 分離）
r_mm = 500.0
theta_deg = 10.0
beam_tan = CurvedMirrorTan(r_mm, theta_deg).apply(beam)  # タンジェンシャル面
beam_sag = CurvedMirrorSag(r_mm, theta_deg).apply(beam)  # サジタル面

print(f'タン: R={beam_tan.R_mm:.1f} mm')
print(f'サグ: R={beam_sag.R_mm:.1f} mm')
```
