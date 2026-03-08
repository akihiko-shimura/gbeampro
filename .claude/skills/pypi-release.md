# Skill: pypi-release

バージョンを上げて PyPI に publish し、GitHub Release を作成する。

## 手順

1. 引数からバージョン種別を判断する（`patch` / `minor` / `major`、デフォルト: `minor`）
2. 現在のバージョンを `gbeampro/__init__.py` から読み取る
3. 新バージョンを計算してユーザーに確認する
4. `gbeampro/__init__.py` の `__version__` を更新する
5. `git add gbeampro/__init__.py` → `git commit -m "Bump version to X.Y.Z"` → `git push`
6. `rm -rf dist && .venv/bin/python -m build` でビルドする
7. `.venv/bin/python -m twine upload dist/*` で PyPI に publish する
8. `gh release create vX.Y.Z dist/*` で GitHub Release を作成する
   - リリースノートには今回の主な変更点を箇条書きで記載する

## 引数

- `patch`: X.Y.Z → X.Y.(Z+1)
- `minor`: X.Y.Z → X.(Y+1).0（デフォルト）
- `major`: X.Y.Z → (X+1).0.0
- バージョン番号を直接指定することも可（例: `2.3.0`）

## 注意

- twine や build が `.venv` にインストールされていない場合はエラーを報告する
- PyPI の認証情報（`~/.pypirc` または環境変数）が必要
