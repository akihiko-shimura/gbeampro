# Project Instructions

## プロジェクト概要

ABCD 行列法（q パラメータ形式）を用いたガウシアンレーザービーム伝播シミュレーション・最適化ライブラリ。
Zemax OpticStudio の Merit Function を参考にした非点収差ビーム整形最適化機能を含む。

## プロジェクト構造

```
gbeampro/       # package source codes
examples/       # Jupyter notebooks of usage examples
docs/           # API reference markdown
.claude/skills/ # project skills
```

## 開発環境

```bash
# 仮想環境
.venv/

# 開発インストール
.venv/bin/pip install -e ".[optimize]"
```

## 言語規約

- `README.md` → 英語
- `docs/*.md` → 日本語

## 利用可能なスキル

| スキル | 用途 |
|--------|------|
| `/nb-exec-push` | notebook 実行 → commit → push |
| `/pypi-release` | バージョンアップ → PyPI → GitHub Release |
| `/api-doc` | Python モジュールから日本語 API reference 生成 |
