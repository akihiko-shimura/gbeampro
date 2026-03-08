# Skill: api-doc

Python モジュールのソースを読んで、日本語の API reference markdown を `docs/` に生成する。

## 手順

1. 引数からモジュール名を特定する（例: `beam`, `elements`, `optimize`）
2. 対応するソースファイル `gbeampro/<module>.py` を読む
3. 既存の `docs/optimize_api.md` のスタイルに倣って `docs/<module>_api.md` を生成する
   - 冒頭にネスト形式の Table of Contents
   - クラス・関数ごとにシグネチャ、パラメータ表、説明
   - コードブロック内コメントは日本語
   - 使用例を含む
4. `docs/` に既存ファイルがあれば上書き前にユーザーに確認する
5. `git add docs/<module>_api.md` → `git commit` → `git push`

## 引数

- モジュール名（単体）: `beam`, `elements`, `system`, `analysis`, `plot`, `optimize`
- 複数指定可: `beam elements` のようにスペース区切り
- 引数なし: 全モジュール分を生成する

## 既存スタイルの参照先

`docs/optimize_api.md` を参照して形式を統一する。

## コミットメッセージ形式

```
Add/Update API reference doc for <module>
```
