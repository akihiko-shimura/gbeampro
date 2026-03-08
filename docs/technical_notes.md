# Technical Notes

## Table of Contents

- [`waist_operands` — curvature weight の注意事項](#waist_operands--curvature-weight-の注意事項)
- [`optimize_astigmatic` — 実行時間の目安](#optimize_astigmatic--実行時間の目安)
- [`TOL_dz_over_zR` — ウェスト位置許容誤差の物理的意味](#tol_dz_over_zr--ウェスト位置許容誤差の物理的意味)
- [`f_abs_bounds` — 非連続探索空間の扱い](#f_abs_bounds--非連続探索空間の扱い)
- [非点収差系の最適化 — ビームパラメータ選択の指針](#非点収差系の最適化--ビームパラメータ選択の指針)

---

## `waist_operands` — curvature weight の注意事項

`waist_tol_mm` による curvature weight 自動計算式:

```
cw = size_weight × (w_mm / (tol_mm / z_mm))²
```

WY が大きく z が大きい場合に `cw` が巨大になり、最適化が失敗しやすい。
例: WY=0.700mm, z=500mm, `waist_tol_y_mm=20` → `cw = (0.700 / (20/500))² = 306`

**推奨**: 不明な場合は `curvature_weight_x=1.0, curvature_weight_y=1.0` を明示的に指定する。

---

## `optimize_astigmatic` — 実行時間の目安

Differential Evolution (`algorithm='de'`) の実行時間:

| `maxiter` | `popsize` | レンズ枚数 | 目安 |
|-----------|-----------|-----------|------|
| 1000 | 15 | 3 | 15–30 s |
| 2000 | 20 | 3 | 30–60 s |

`algorithm='hammer'` は DE + TRF 精密化のため上記の 1.5–2 倍程度。

---

## `TOL_dz_over_zR` — ウェスト位置許容誤差の物理的意味

PASS 条件: `|Δz| / z_R ≤ TOL_dz_over_zR`

- `z_R = π w₀² / λ` はビームウェストの Rayleigh 長
- `TOL_dz_over_zR = 0.2` は Rayleigh 長の 20% 以内にウェストが来ることを要求
- ビームが太い（w₀ 大）ほど z_R が大きくなるため、絶対量としては緩い制約になりうる
- 逆に、強収束（w₀ 小）では z_R が短く、制約が厳しくなる

---

## `f_abs_bounds` — 非連続探索空間の扱い

`f_abs_bounds=(f_lo, f_hi)` を指定すると許容範囲は `[-f_hi, -f_lo] ∪ [f_lo, f_hi]`。

DE はこの非連続空間を直接扱えないため、`_residuals` 内でペナルティを付与する実装になっている:

```python
def _fabs_ok(params):
    fs = [_snap_f(f) for f in params[1::2]]
    return all(abs(f) >= f_abs_lo for f in fs)

def _residuals(params):
    if not _sep_ok(params) or not _fabs_ok(params):
        return np.full(len(operands), 1e3)  # penalty
    ...
```

`f_step_mm > 0` の場合、スナップ後の値で判定する必要があるため `_snap_f` を先に適用すること。

---

## 非点収差系の最適化 — ビームパラメータ選択の指針

x/y 軸のターゲット比 (WX : WY) と球面・シリンドリカルレンズの組み合わせ:

- WX = WY（円形収束）→ 球面レンズのみで原理上解ける
- WX ≠ WY（非点収差）→ `spherical + cyl_x + cyl_y` の 3 枚構成が基本
- 収束比 `w₀ / W0` が大きいほど（強収束）短焦点レンズが必要になり、`f_abs_bounds` の下限制約が効いてくる
