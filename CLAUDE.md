# Project Instructions

## プロジェクト構造

```
gbeampro/       # パッケージ本体
examples/       # Jupyter notebooks（プロジェクトルートに配置）
docs/           # API reference markdown
.claude/skills/ # プロジェクト固有のスキル
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

## リリース手順

1. `gbeampro/__init__.py` の `__version__` を更新
2. `git commit` → `git push`
3. `rm -rf dist && .venv/bin/python -m build`
4. `.venv/bin/python -m twine upload dist/*`
5. `gh release create vX.Y.Z dist/*`

→ `/pypi-release` スキルで自動化可能

## 利用可能なスキル

| スキル | 用途 |
|--------|------|
| `/nb-exec-push` | notebook 実行 → commit → push |
| `/pypi-release` | バージョンアップ → PyPI → GitHub Release |
| `/api-doc` | Python モジュールから日本語 API reference 生成 |
