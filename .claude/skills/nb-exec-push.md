# Skill: nb-exec-push

Jupyter notebook を実行して出力を埋め込み、結果を確認してから commit & push する。

## 手順

1. 対象の notebook のパスを特定する（引数があればそれを使う。なければ最近編集された notebook を探す）
2. `CLAUDE.md` の手順に従い、`.venv/bin/jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=600 <path>` を実行する
3. 実行結果（主要セルの出力）を確認してユーザーに報告する
4. エラーがあればユーザーに報告して停止する
5. 問題なければ `git add <path>` → `git commit` → `git push` を実行する

## 引数

- 引数あり: 指定された notebook パスを使う
- 引数なし: 最近変更された `.ipynb` ファイルを対象にする

## コミットメッセージ形式

```
Re-execute <notebook_name>: <変更内容の一言説明>
```
