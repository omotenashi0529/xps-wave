# xps-wave

XPS（X線光電子分光）測定データの波形分離（ピークデコンボリューション）、ピーク位置（結合エネルギー）・強度の算出、元素・化学状態推定までを行う Python パッケージ。Ti (チタン) 2p 領域の高精度解析と、ワイドスペクトル（サーベイスキャン）の自動解析の両方に対応する。

## セットアップ (uv)

```bash
uv sync          # 依存関係をインストール (numpy, scipy, lmfit, pandas, matplotlib, vamas, requests)
uv run pytest    # ユニットテスト (11件) を実行
```

Python 3.12 以上、`pyproject.toml` に依存関係を固定済み。仮想環境は `.venv/` に作成される。

## クイックスタート

```python
from xps_wave.io import read_vamas, read_two_column
from xps_wave.pipeline import XPSAnalyzer
from xps_wave.plotting import plot_analysis

# CasaXPS / Kratos / PHI などが出力する VAMAS (.vms) ファイル
spectrum = read_vamas("data/raw/mixed_ti_2p.vms")

# もしくは 2 列 (エネルギー, 強度) のテキスト/CSV
# spectrum = read_two_column("my_region.csv", delimiter=",", skip_header=1)

analyzer = XPSAnalyzer(spectrum)

# Ti 2p 領域: 文献値で制約したスピン軌道二重線モデルで解析
result = analyzer.analyze_ti_2p(low=450, high=470)
print(result.table())   # BE(eV), 強度(Y軸), 面積, FWHM, 推定元素/化学状態, 文献値との誤差

plot_analysis(result, title="Ti 2p").figure.savefig("ti2p.png")
```

`result.table()` の出力例（実測データでの実行結果）:

```
                      label   BE_eV  intensity     area  fwhm_eV element orbital                 state  ref_BE_eV  delta_eV confidence
             Ti metal 2p3/2 453.720     5358.9   5436.9    0.953      Ti   2p3/2              Ti metal     453.80    -0.080       high
                  TiO 2p3/2 454.600     2822.9   7378.6    2.456      Ti   2p3/2                   TiO     455.40    -0.800        low
                Ti2O3 2p3/2 457.494     2561.0   5083.2    1.865      Ti   2p3/2                 Ti2O3     457.20     0.294       high
TiO2 (rutile/anatase) 2p3/2 458.606    14847.5  17800.7    1.126      Ti   2p3/2 TiO2 (rutile/anatase)     458.60     0.006       high
             Ti metal 2p1/2 459.620     2250.5   2718.4    1.135      Ti   2p1/2              Ti metal     459.85    -0.230       high
                  ...
```

## アーキテクチャ

```
src/xps_wave/
  spectrum.py    Spectrum データクラス、KE<->BE 変換
  io.py          読み込み: 汎用2列テキスト/CSV、VAMAS(.vms)。CasaXPS の校正情報(Calib M=.. A=..)を自動検出・適用
  background.py  バックグラウンド除去: Shirley（逐次法）、線形、ローリング最小値（サーベイ用）
  peakfit.py     波形分離 (lmfit PseudoVoigt): 汎用マルチピーク fit_peaks / Ti 2p 専用スピン軌道二重線 fit_ti_2p_doublets
  reference.py   結合エネルギー参照データベース（Ti 2p 高精度版 + 30元素程度の広域テーブル、出典明記）
  identify.py    フィッティング結果を参照データベースと照合し元素・化学状態を推定
  survey.py      ワイドスペクトル全域のピーク自動探索・個別フィット・同定
  pipeline.py    XPSAnalyzer: 上記を束ねるエンドツーエンドの解析クラス
  plotting.py    生データ・バックグラウンド・フィット結果・分離ピークの可視化
scripts/
  fetch_reference_data.py   実測 Ti 2p 参照スペクトルをインターネットから取得
  evaluate_accuracy.py      実測データに対する精度評価（本パイプライン vs 専門家フィット）
  analyze_wide_spectrum.py  ワイドスペクトル解析のデモ（合成スペクトル、下記参照）
tests/           background/peakfit/identify のユニットテスト
data/raw/        取得した実測データ
results/         精度評価・デモの出力（CSV, PNG）
```

## 解析アルゴリズムの要点

- **バックグラウンド除去**: 標準的な Shirley 法（逐次計算、`background.shirley_background`）。線形法、サーベイ用のローリング最小値法も選択可能。
- **波形分離**: lmfit の `PseudoVoigtModel`（Gaussian-Lorentzian混合）を各ピークに割り当てて同時フィット。
  - Ti 2p では化学状態ごとに 2p3/2 / 2p1/2 のスピン軌道二重線を構成し、
    - 面積比を 4:2 (縮重度 2:1) に固定（自由パラメータにしない — 原子物理で決まる値のため）
    - 分裂幅を文献値 ± 0.15 eV の範囲でのみ微調整可能
    - 幅(sigma)を文献 FWHM の 0.5〜1.8 倍に制限（重なりの強い複数化学状態が縮退して1本に潰れる問題を回避）
    という物理的制約を組み込むことで、汎用フィットより大幅に安定した結果を得ている。
- **元素推定**: フィットしたピーク位置(BE)を `reference.py` の参照データベースと照合し、最も近い候補（許容誤差内、既定 ±1.0 eV）を第一候補として提示。Ti 2p 解析では二重線モデル自体が構造的に化学状態を確定させているため、参照値との誤差(delta_eV)は「クロスチェック」として使う。

## データ出所・精度評価

**インターネットから取得した実測データ**: [XPS Reference Pages: Titanium](http://www.xpsfitting.com/2008/09/titanium.html)（M.C. Biesinger, Surface Science Western, Western University）で公開されている "Mixed Titanium Sample, Ti 2p" の VAMAS (.vms) ファイル。金属Ti・TiO・Ti2O3・TiO2 を混合した実試料の実測スペクトルで、CasaXPS による専門家フィット結果（8ピークの結合エネルギー・面積・幅）がファイル内にテキストとして埋め込まれている。`scripts/fetch_reference_data.py` が実行時にダウンロードする（リポジトリには同梱していない）。

> 取得時の注意点: このファイルの生スペクトル軸は未校正（as-measured）の運動エネルギー軸で、校正情報（`Calib M=455.59 A=458.6 BE ADD`）はコメントとして別途埋め込まれていた。これに気づかず素朴に KE→BE 変換すると主成分の帰属が Ti(II) と Ti(IV) で入れ替わってしまう（`io.py` の `read_vamas` はこの校正シフトを自動検出・適用する）。実際の XPS 解析でも帯電補正（charge referencing）は必須であることを示す実例。

**精度評価** (`uv run python scripts/evaluate_accuracy.py`): このファイルに埋め込まれた専門家フィット（CasaXPS）を正解値として、本パッケージ自身のフィッティング結果と比較。

| 指標 | 値 |
|---|---|
| 全8ピーク平均誤差 (MAE) | 0.30 eV |
| RMSE | 0.46 eV |
| 最大誤差 | 0.95 eV (Ti 2p1/2 Ti(II)) |
| 主成分 (Ti metal, TiO2) の MAE | **0.05 eV** |
| 少数・重なり成分 (TiO, Ti2O3) の MAE | 0.55 eV |

主成分（金属Ti、TiO2）は 0.05 eV というXPSの実用精度として十分なレベルで再現できる。少数成分の Ti(II)/Ti(III) はスペクトル的に強く重なり合っており（そもそも Biesinger らの原論文がこの問題を解くために書かれている）、誤差が大きくなるのは物理的に妥当な結果であり、本パイプラインの限界として明記する。結果は `results/ti2p_accuracy.png`（可視化）と `results/ti2p_accuracy_report.csv` に保存される。

## ワイドスペクトル対応

`xps_wave.survey.analyze_wide_spectrum()` で全域スキャンから自動でピークを検出し、`reference.SURVEY_REFERENCE_TABLE`（約30元素の主要ピーク）と照合して元素を推定する。Ti 領域に絞りたい場合は `XPSAnalyzer.analyze_ti_2p()` を使う方が精度が高い（上記の物理制約付きモデルのため）。

`uv run python scripts/analyze_wide_spectrum.py` で動作を確認できる。ただし、このデモは**合成（synthetic）スペクトル**を使用している — ライセンスが明確でオープンに再配布可能な実測ワイドスペクトルは、本プロジェクト構築時の調査では見つからなかった（詳細は同スクリプトのdocstring参照）。合成スペクトルは文献値の結合エネルギー・現実的な幅・Poisson統計ノイズ・段差状バックグラウンドを組み込んで物理的に妥当な形にしているが、これは「精度評価」ではなく「コードパスの動作確認」を目的としたものである点に注意。実際のワイドスペクトルを使う場合は `read_two_column()` または `read_vamas()` で読み込んだ実データをそのまま渡せる。

## 既知の限界

- `reference.SURVEY_REFERENCE_TABLE` の値は NIST XPS Database / Moulder Handbook 等でよく引用される代表的な文献値であり、定量的な厳密性が必要な場合は [NIST XPS Database](https://srdata.nist.gov/xps/) を直接参照すること。
- 帯電補正（charge referencing、例: C1s adventitious carbon = 284.8 eV への補正）は自動化していない。ユーザー側で校正済みのBEスケールにしてから読み込むか、`Spectrum.energy` を直接シフトする必要がある（VAMASファイルにCasaXPS校正情報が埋め込まれている場合は自動適用される）。
- ワイドスペクトルの単一ピークフィットは、スピン軌道二重線などの重なりを積極的にはモデル化しない（Ti は専用の `analyze_ti_2p()` を推奨）。

## 推奨する Skill / MCP

このタスクに現時点でジャストフィットする既存の Skill/MCPは無かったため、以下を提案する。

1. **カスタム Skill "xps-analysis" の作成（推奨）**: 本リポジトリの `scripts/evaluate_accuracy.py` 相当のワークフロー（データ読み込み→背景除去→フィット→同定→レポート）を `SKILL.md` として切り出せば、次回以降 "このXPSデータを解析して" と頼むだけで一連の処理を再現できる。`find-skills` Skillで類似の既存Skillが公開されていないか確認可能。
2. **Jupyter/ノートブック実行系 MCP**: `mcp__ide__executeCode` は既に利用可能（このセッションでも使用したPython実行環境の元）。装置から出力されたVAMASファイルを対話的に読み込んで即座にプロットしながら調整したい場合に有用。
3. **文献・データベース検索用 MCP（あれば）**: NIST XPS Database や CCDC/PubChem 等、材料科学系のデータベースに直接クエリできるMCPサーバーがあれば、`reference.py` の値をハードコードでなく動的に検証・拡充できる。現時点で汎用的にAPI化されたNIST XPS DatabaseのMCP実装は確認できなかった（Web検索インターフェースのみ）。
4. 既存の汎用 Skill としては **dataviz**（グラフ・可視化のデザイン指針）が `plotting.py` の見た目を改善する際に役立つ。

## テスト

```bash
uv run pytest -q
```

11件のユニットテスト（Shirley/線形バックグラウンドの境界条件、単一・重複ガウシアンのピーク分離、スピン軌道面積比の物理制約、元素同定の信頼度判定）が独立して実行できる。実測データを使った統合的な精度評価は `scripts/evaluate_accuracy.py`（インターネット接続が必要）で行う。
