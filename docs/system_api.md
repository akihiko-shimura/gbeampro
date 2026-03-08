# `gbeampro.system` API Reference

## Table of Contents

- [Overview](#overview)
- [クラス](#クラス)
  - [`OpticalSystem`](#opticalsystem)
    - [`add`](#add)
    - [`trace`](#trace)
    - [`__str__`](#__str__)
    - [`summary`](#summary)
- [使用例](#使用例)
  - [光学系の構築とビーム追跡](#光学系の構築とビーム追跡)
  - [素子レイアウトの表示](#素子レイアウトの表示)
  - [ウェスト情報を含むサマリーの表示](#ウェスト情報を含むサマリーの表示)

---

## Overview

複数の光学素子を連結し，ガウシアンビームを追跡するクラス `OpticalSystem` を提供する。

- **Fluent API** により素子を連鎖的に追加できる
- `trace` は `Propagation` 素子内の中間点も含む全軌跡を返す
- `summary` は各素子境界でのビーム状態とウェスト情報を人間が読める形式で出力する

---

## クラス

### `OpticalSystem`

光学素子のシーケンスを管理し，ビーム追跡を行うクラス。

```python
from gbeampro.system import OpticalSystem

system = OpticalSystem()
```

---

### `add`

素子を光学系の末尾に追加する。`self` を返すため，メソッドチェーンで記述できる。

```python
def add(self, element: Element) -> OpticalSystem
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `element` | 追加する `Element` インスタンス |

**戻り値**: `self`（fluent API）

```python
system = (OpticalSystem()
    .add(Propagation(100))
    .add(ThinLens(f_mm=150))
    .add(Propagation(200))
)
```

---

### `trace`

入力ビームをシステム全体に伝搬させ，全軌跡を `GaussBeam` のリストで返す。

```python
def trace(self, beam: GaussBeam, dz: float = 0.01) -> list[GaussBeam]
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `beam` | 入力ビーム（初期状態） |
| `dz` | `Propagation` 素子内で中間状態を記録する刻み幅 (mm) |

**戻り値**: `list[GaussBeam]`

- リストの各要素は `z_mm`, `w_mm`, `R_mm` を持つ `GaussBeam`
- 先頭要素は入力ビームそのもの（`beam`）
- `Propagation` 素子は `dz` 刻みで中間状態を含む（caustic プロットに使用）
- `ThinLens`, `Interface` など厚みのない素子は1点のみ追加

**中間点の含め方**

`Propagation(d_mm=100)` に `dz=1.0` を指定すると，`z=0, 1, 2, ..., 100` の計 101 点が軌跡に含まれる。`dz` を小さくするほどプロットが滑らかになるが，計算コストが増す。

---

### `__str__`

素子レイアウトを表形式の文字列で返す。`print(system)` で確認できる。

```python
def __str__(self) -> str
```

**出力例**

```
#   Type            z_mm     param
--  --------------  -------  --------
 0  Propagation     100.0    d=100.0
 1  ThinLens        100.0    f=150.0
 2  Propagation     300.0    d=200.0
```

---

### `summary`

各素子境界でのビーム状態（z, w, R, θ）とウェスト位置・半径の一覧を文字列で返す。

```python
def summary(self, beam: GaussBeam, dz: float = 0.01) -> str
```

**パラメータ**

| 引数 | 説明 |
|------|------|
| `beam` | 入力ビーム |
| `dz` | `trace` に渡す刻み幅 (mm) |

**出力に含まれる情報**

1. **素子境界でのビーム状態テーブル**: 各素子の前後における z, w (mm), R (mm), θ (mrad)
2. **ウェスト情報**: 検出されたウェスト位置 z_mm と半径 w_mm の一覧

**出力例**

```
Beam state at each element boundary:
  z_mm     w_mm      R_mm      theta_mrad
  -------  --------  --------  ----------
    0.000    1.0000       inf      0.3396
  100.000    1.0341   3030.30      0.3396
  100.000    1.0341   -170.45      0.3396
  300.000    0.1823       inf      0.3396

Waists found:
  z_mm=300.00  w_mm=0.1823
```

---

## 使用例

### 光学系の構築とビーム追跡

```python
from gbeampro import GaussBeam
from gbeampro.elements import Propagation, ThinLens
from gbeampro.system import OpticalSystem

# 入力ビーム: λ=1064 nm, ウェスト半径 1.0 mm
beam = GaussBeam.from_waist(wl_um=1.064, w0_mm=1.0)

# 光学系: 100 mm 伝搬 → f=150 mm レンズ → 200 mm 伝搬
system = (OpticalSystem()
    .add(Propagation(100))
    .add(ThinLens(f_mm=150))
    .add(Propagation(200))
)

# ビーム軌跡を取得（dz=1.0 mm 刻み）
trajectory = system.trace(beam, dz=1.0)

# 各点でのビーム半径を確認
for b in trajectory[::50]:  # 50 点おきにサンプリング
    print(f'z={b.z_mm:.1f} mm,  w={b.w_mm:.4f} mm')
```

### 素子レイアウトの表示

```python
from gbeampro.elements import Propagation, ThinLens, Interface
from gbeampro.system import OpticalSystem

system = (OpticalSystem()
    .add(Propagation(50))
    .add(Interface(n1=1.0, n2=1.5))   # ガラス入射界面
    .add(Propagation(30))
    .add(Interface(n1=1.5, n2=1.0))   # ガラス出射界面
    .add(Propagation(100))
    .add(ThinLens(f_mm=200))
    .add(Propagation(300))
)

print(system)  # 素子レイアウト表を出力
```

### ウェスト情報を含むサマリーの表示

```python
from gbeampro import GaussBeam
from gbeampro.elements import Propagation, ThinLens
from gbeampro.system import OpticalSystem

beam = GaussBeam.from_waist(wl_um=0.8, w0_mm=0.5)

system = (OpticalSystem()
    .add(Propagation(200))
    .add(ThinLens(f_mm=100))
    .add(Propagation(200))
)

# 素子境界でのビーム状態 + ウェスト情報を一括表示
print(system.summary(beam, dz=0.5))
```
