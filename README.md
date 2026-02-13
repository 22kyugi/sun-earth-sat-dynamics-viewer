# sun-earth-sat-dynamics-viewer

太陽中心（ヘリオセントリック）視点で地球の公転と衛星の地球周回を可視化する MATLAB スクリプトです。
UI（スライダー）で月 / 日 / 時刻を操作し、任意時刻の太陽方向・昼夜境界（ターミネータ）・衛星位置・軌道面を表示できます。

English (short): A MATLAB visualization tool for Sun-centered Earth orbit and Earth-centered satellite motion with UI sliders (month/day/hour). It renders satellite CAD models (STL), supports attitude modes (nadir/sun/velocity/keyframes), and provides playback controls and a satellite zoom inset with an attitude triad.

---

## Features

- Sun-centered view (left)
  - 地球の公転軌道（1年分）
  - 地球の移動（テクスチャ付き地球、昼夜陰影、ターミネータ）
  - 春分・夏至・秋分・冬至のマーカー表示（近似日付）
  - 太陽の橙色マーカー表示（原点）

- Earth-centered view (right)
  - 地球テクスチャに地球自転（GMST）を反映（テクスチャが回転）
  - 衛星の現在位置（STLモデル表示 / フォールバック点）
  - 現在時刻の瞬時軌道（軌道線）のみを表示
  - 軌道面（半透明ディスク） + 法線ベクトル矢印
  - 簡易 地球影（食）判定（円柱影近似）

- Satellite Zoom Inset (right-top)
  - 右図右上に衛星のみの拡大図（インセット）を表示（地球は表示しない）
  - 拡大図に衛星中心からのXYZ矢印（トライアド）を表示（衛星姿勢に追従）
  - 矢印長は拡大図のズーム幅に応じて自動調整（見切れにくい）
  - インセットはサイズ調整可能（軌道線への被りを軽減するため縮小・右上詰め設定）

- Attitude (姿勢)
  - 姿勢モード：fixed / keyframes / nadir / sun / velocity をサポート
  - nadir/sun/velocity では軸割当（xAxis/zAxis）+ roll（rollDeg/rollAxis）を指定可能
  - オプションで区間ごとの姿勢モード切替（スケジュール）に対応

- UI
  - スライダー：月 → 日 → 時刻（時刻は分刻み量子化に対応）
  - 右側に Play / Pause を配置し、12h/24h（半日/1日）ステップで自動進行可能
  - 進行速度（タイマー周期）を調整可能（Speed[s]）

- STL loading
  - stlread が使える環境ではそれを利用
  - stlread が無い場合でも内蔵STLリーダ（ASCII/Binary）で読み込みを試行
  - 読み込み失敗時はフォールバック点表示

---

## Demo / Screenshot

![スクリーンショット](docs/screenshot.png)

図の説明（左：太陽中心 / 右：地球中心）
- 左図（太陽中心）：太陽を原点として、地球の公転軌道（1年分）と現在位置（青点）を表示します。地球にはテクスチャと昼夜陰影（ターミネータ）が描かれ、季節変化が直感的に分かります。
- 右図（地球中心）：衛星の現在位置（STLモデル）と瞬時軌道（軌道線）、軌道面（半透明ディスク）を表示します。地球テクスチャはGMSTに基づいて回転します。
- 右上拡大図：衛星のみをズーム表示し、衛星中心からXYZ矢印（トライアド）を表示します。矢印は衛星姿勢に追従し、長さはズーム幅に応じて自動調整されます。
- 下部UI：スライダーで日時を指定、Play/Pauseで12h/24hステップの自動進行が可能です。

---

## Usage

1. visualizeSunEarthSat_UI_3sliders.m を MATLAB のカレントフォルダへ配置
2. （任意）satellite.stl を同じフォルダへ配置（または params.satModelFile をフルパス指定）
3. 実行：

```matlab
visualizeSunEarthSat_UI_3sliders
```
---
## Parameters
主要なパラメータは visualizeSunEarthSat_UI_3sliders.m 冒頭（User Parameters）で変更できます。
詳細は PARAMETERS.md を参照してください。

---
## Troubleshooting
STLが点表示になる（モデルが表示されない）

- STLファイルが読み込めていない可能性があります。
- satellite.stl がカレントフォルダにあるか確認してください。
- params.satModelFile をフルパス指定してみてください。
```matlab
params.satModelFile = "C:\path\to\satellite.stl";
```
---
## License

This project is licensed under the MIT License.