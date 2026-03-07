# Project Instructions

## Jupyter Notebook の push 手順

notebook を push する前に必ず以下の手順を踏む：

1. VS Code でファイルを保存済みであることを確認する
2. `jupyter nbconvert --to notebook --execute --inplace` で実行結果を埋め込む
3. 出力が正しく含まれていることを確認してから `git add` → `commit` → `push`

```bash
.venv/bin/jupyter nbconvert --to notebook --execute --inplace <notebook_path>
```

実行に失敗した場合はエラーを確認し、ユーザーに報告してから push しない。
