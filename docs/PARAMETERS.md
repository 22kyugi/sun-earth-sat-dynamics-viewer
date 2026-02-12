# PARAMETERS.md — パラメータ解説

このドキュメントは `visualizeSunEarthSat_UI_3sliders.m` の `User Parameters` で調整可能な設定の説明です。  
特に今回追加・調整した **インセット（拡大図）** と **XYZ矢印（トライアド）**、**時間刻み**のパラメータを中心にまとめます。

---

## 1. 時間・UI関連

### `params.yearOffsetN`
- **意味**：基準年を `2000 + yearOffsetN` で指定
- **例**：`26` → 2026年

### `params.timeStepMin`
- **意味**：Hourスライダーの値を「何分刻みで量子化するか」
- **単位**：分
- **例**
  - `1` → 1分刻み（最も細かい）
  - `5` → 5分刻み
  - `10` → 10分刻み
- **挙動**：スライダー値を `round(value/step)*step` の形で丸めて時刻に反映するため、表示時刻は必ず指定刻みに揃います。

---

## 2. インセット（拡大図）関連

### `params.useInsetZoom`
- **意味**：右上の拡大図（インセット）を表示するか
- **true/false**

### `params.insetSizeFactor`
- **意味**：右図（ax2）のタイルに対するインセットのサイズ比率
- **例**
  - `0.30` → 右図の幅/高さの30%
  - `0.15` → 約半分（軌道線への被りが減る）
- **推奨**：被りが気になる場合 `0.12 ～ 0.18`

### `params.insetMargin`
- **意味**：インセットを右上に配置する際の余白（normalized座標）
- **例**
  - `0.010` → 余白やや大きめ
  - `0.002` → 右上へ詰める

### `params.insetHalfWidth_km`
- **意味**：拡大図の表示範囲（衛星中心 ±dz）
- **単位**：km
- **例**：`20000` → 衛星中心からX/Y/Zそれぞれ ±20000 km の範囲を表示
- **調整指針**
  - 衛星を大きく表示している場合は、枠内に収まるよう `dz` を広げる
  - 衛星をより拡大したい場合は `dz` を小さくする

> 拡大図は「衛星のみ」表示のため、衛星が地球から離れると地球は表示されません（意図通り）。

---

## 3. 衛星STLモデル表示関連

### `params.useSatModel`
- **意味**：STLモデルを使用するか
- **true/false**
- falseの場合は点表示（フォールバック）

### `params.satModelFile`
- **意味**：STLファイル名またはパス
- **例**
  - `"satellite.stl"`（カレントフォルダに配置）
  - `"C:\path\to\satellite.stl"`

### `params.satModelUnit`
- **意味**：STLモデルの単位
- **値**：`"mm"`, `"m"`, `"km"`
- **注意**：STLは単位情報を持たないため、ここで正しく指定する必要があります。

### `params.satModelDesiredSize_km_main`
- **意味**：右図メイン（Earth-centered）での衛星モデルの「表示サイズ目標」
- **単位**：km  
- **注意**：大きくしすぎると地球や軌道線を覆います。

### `params.satModelDesiredSize_km_inset`
- **意味**：拡大図での衛星モデルの「表示サイズ目標」
- **単位**：km  
- **用途**：拡大図は衛星確認用なので大きめに設定してOK

### `params.fallbackMarkerColor`
- **意味**：STLが読めない場合に表示する点（フォールバック）の色
- **例**：`[0.1 0.1 0.1]`（黒系）

---

## 4. 拡大図のXYZ矢印（トライアド）関連

### `params.showInsetTriad`
- **意味**：拡大図にXYZ矢印を表示するか
- **true/false**

### `params.triadLineWidth`
- **意味**：矢印の線幅

### `params.triadLength_km`
- **意味**：矢印長の上限値（km）
- **実際の挙動**：自動調整後に `min( triadLength_km, triadAutoFrac*dz )` を適用  
  → `triadLength_km` は「最大値」として働きます

### `params.triadAutoFrac`
- **意味**：インセット表示範囲 `dz` に対する矢印長の割合
- **式**：`Lauto = triadAutoFrac * dz`
- **例**
  - `0.35` → 枠の35%程度の矢印長（見切れにくい）

### `params.showTriadLabels`
- **意味**：矢印先端に X/Y/Z のラベルを表示するか
- **true/false**
- **挙動**：ラベル位置も `rSat + axisVector` で更新し、姿勢に追従します。

---

## 5. 推奨設定例

### 被りを減らす（インセット小さめ＋右上詰め）
```matlab
params.insetSizeFactor = 0.15;
params.insetMargin     = 0.002;