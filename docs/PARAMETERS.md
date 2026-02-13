# PARAMETERS.md - パラメータ解説

このドキュメントは visualizeSunEarthSat_UI_3sliders.m の User Parameters で調整できる設定の説明です。
本プログラムは以下の機能を含みます。

## 目次

- [PARAMETERS.md - パラメータ解説](#parametersmd---パラメータ解説)
  - [目次](#目次)
  - [01. 年・時刻・UI関連](#01-年時刻ui関連)
  - [02. 再生（Play/Pause）関連](#02-再生playpause関連)
  - [03. 軌道要素（衛星軌道）関連](#03-軌道要素衛星軌道関連)
  - [04. 左図（太陽中心表示）と右図（地球中心表示）関連](#04-左図太陽中心表示と右図地球中心表示関連)
  - [05. 昼夜陰影（地球テクスチャ/照明）関連](#05-昼夜陰影地球テクスチャ照明関連)
  - [06. 地球自転（GMST）関連](#06-地球自転gmst関連)
  - [07. 右上インセット（衛星のみの拡大表示）関連](#07-右上インセット衛星のみの拡大表示関連)
  - [08. インセット上のXYZ矢印（トライアド）関連](#08-インセット上のxyz矢印トライアド関連)
  - [09. 姿勢（Attitude）モード関連](#09-姿勢attitudeモード関連)
    - [09.1 fixed モード](#091-fixed-モード)
    - [09.2 keyframes モード](#092-keyframes-モード)
    - [09.3 nadir / sun / velocity モード（指向＋ロール）](#093-nadir--sun--velocity-モード指向ロール)
  - [10. 区間スケジュール（姿勢モード切替）オプション](#10-区間スケジュール姿勢モード切替オプション)
  - [11. STLモデル（衛星形状）表示関連](#11-stlモデル衛星形状表示関連)
  - [12. 推奨設定例](#12-推奨設定例)
  - [13. トラブルシューティング（よくある症状）](#13-トラブルシューティングよくある症状)

---

## 01. 年・時刻・UI関連

params.yearOffsetN
- 意味：基準年を 2000 + yearOffsetN で指定する
- 例：26 -> 2026年

params.baseYear
- 意味：実際に使う年（baseYear = 2000 + yearOffsetN）
- 注意：yearOffsetN を変えると自動的に変わる想定

params.timeStepMin
- 意味：Hourスライダーの値を「何分刻み」で量子化する
- 単位：分
- 例：1 / 5 / 10
- 補足：1分刻みにすると操作は細かくなるが、更新頻度が上がるため環境によっては重くなることがある

---

## 02. 再生（Play/Pause）関連

この機能は UI 右側の Play/Pause と Step, Speed で制御されます。
一部の値は内部状態変数として持ちます（User Parameters ではなく関数内の状態）。

playStepDays（内部状態）
- 意味：再生時に進める時間ステップ
- 値：0.5（12時間） / 1.0（24時間）
- UI：Step の 12h/24h で切替

playPeriodSec（内部状態）
- 意味：タイマー周期（何秒ごとに次ステップへ進めるか）
- 単位：秒
- UI：Speed[s] の入力値
- 目安：0.2～1.0 秒程度（PC性能に応じて調整）

currentDT（内部状態）
- 意味：現在表示中のUTC時刻（再生の基準時刻）
- 補足：スライダー操作でも更新される（再生と手動が同期）

isAutoUpdating（内部状態）
- 意味：再生がスライダー値を書き換える時の再入防止フラグ
- 目的：スライダー更新のコールバックが再帰的に呼ばれて無限ループになるのを防ぐ

注意
- 図を閉じる時に timer を停止・削除しないと、MATLABセッションに timer が残る場合があるため、CloseRequestFcn で停止・削除する実装が推奨

---

## 03. 軌道要素（衛星軌道）関連

params.orbit.a
- 意味：長半径
- 単位：km
- 例：Re + alt_km

params.orbit.e
- 意味：離心率

params.orbit.i
- 意味：軌道傾斜角
- 単位：deg

params.orbit.RAAN
- 意味：昇交点赤経
- 単位：deg

params.orbit.argp
- 意味：近地点引数
- 単位：deg

params.orbit.M0
- 意味：平均近点角（初期）
- 単位：deg

params.orbit.useJ2
- 意味：J2による歳差（RAAN, argp）の簡易モデルを入れるか
- true/false
- 注意：可視化用の簡易モデルで、厳密な軌道力学モデルではない

---

## 04. 左図（太陽中心表示）と右図（地球中心表示）関連

params.AU_scaleDraw
- 意味：AUスケールの描画倍率（AU距離の表示倍率）
- 通常：1

params.sunDrawRadius_km
- 意味：描画用の太陽半径（見やすさ用に拡大）
- 単位：km

params.earthDrawRadius_km
- 意味：描画用の地球半径（見やすさ用に拡大）
- 単位：km

params.satRelScale
- 意味：太陽中心表示での衛星位置を地球周りで見えるように拡大する倍率（描画専用）
- 注意：物理スケールとは一致しない（見やすさ重視）

---

## 05. 昼夜陰影（地球テクスチャ/照明）関連

params.useTexture
- 意味：地球表面をテクスチャ表示にするか
- true/false
- false の場合は陰影のみの簡易表現

params.nightFactor
- 意味：夜側の明るさ（0に近いほど暗い）
- 範囲：0～1
- 例：0.12

params.gamma
- 意味：照明のガンマ調整（陰影の強調）
- 例：0.75

---

## 06. 地球自転（GMST）関連

右図（地球中心）では地球テクスチャをGMSTに基づいて回転させます。
調整用パラメータは通常不要で、dtUTCから内部計算します。

---

## 07. 右上インセット（衛星のみの拡大表示）関連

params.useInsetZoom
- 意味：右上インセットを表示するか
- true/false

params.insetSizeFactor
- 意味：右図タイルに対するインセットサイズ比率（幅/高さに掛かる）
- 例：0.15（小さめ、被り軽減）

params.insetMargin
- 意味：インセットを右上に詰める際の余白（normalized座標）
- 例：0.002（かなり詰める）

params.insetHalfWidth_km
- 意味：インセットの表示範囲（衛星中心 ±dz）
- 単位：km
- 値を小さくするとよりズーム、値を大きくすると広域表示

---

## 08. インセット上のXYZ矢印（トライアド）関連

params.showInsetTriad
- 意味：インセットにXYZ矢印を表示するか
- true/false

params.triadAutoFrac
- 意味：インセット表示範囲 dz に対する矢印長の割合
- 例：0.35
- 計算：Lauto = triadAutoFrac * dz

params.triadLength_km
- 意味：矢印長の上限（km）
- 実際の矢印長は min(triadLength_km, Lauto) で決まる
- 補足：小さすぎる場合の下限（例：0.05*dz）を追加して視認性を確保する実装が推奨

params.triadLineWidth
- 意味：矢印線幅

params.showTriadLabels
- 意味：矢印先端に X/Y/Z ラベルを表示するか
- true/false

---

## 09. 姿勢（Attitude）モード関連

姿勢は姿勢行列 Rbody（body -> ECI）として扱い、STLモデルとトライアドに適用します。

params.attMode
- 意味：姿勢モードを選択する
- 値：fixed / keyframes / nadir / sun / velocity

params.satSpinRPM
- 意味：姿勢適用後にボディZ軸周りに回すスピン
- 単位：RPM
- 補足：固定姿勢や指向姿勢の上に追加回転として適用される

### 09.1 fixed モード

params.satEuler_deg
- 意味：固定Euler角（deg）
- 例：[140 330 90]

params.satEulerOrder
- 意味：Euler角の回転順序
- 例：XYZ / ZYX

### 09.2 keyframes モード

params.attKF
- 意味：時刻とEuler角を持つ姿勢スケジュール（table）
- 想定列：tUTC, roll, pitch, yaw
- 注意：datetime と数値を同じ配列に入れず、table で保持する

params.attInterp
- 意味：補間方式
- 例：slerp

params.attEulerOrder
- 意味：keyframes のEuler角回転順序

### 09.3 nadir / sun / velocity モード（指向＋ロール）

params.attOpt
- 意味：指向軸とロールを指定する構造体

attOpt.xAxis
- 意味：ボディX軸を向ける方向
- 値例：vel / sun / nadir / orbitNormal / antivel / antisun など
- 実装が対応する spec のみ有効

attOpt.zAxis
- 意味：ボディZ軸を向ける方向（xAxisと直交化してフレーム構成）
- 値例：nadir / sun / orbitNormal など

attOpt.rollDeg
- 意味：ロール角（deg）

attOpt.rollAxis
- 意味：どのボディ軸周りにロールするか
- 値：x / y / z

指向フレームの考え方
- まず xAxis を厳密に合わせる
- zAxis は xAxis と直交化してできるだけ希望方向に近づける
- yAxis は外積で生成して直交系にする

---

## 10. 区間スケジュール（姿勢モード切替）オプション

params.useAttSchedule
- 意味：区間ごとの姿勢モード切替を有効にする
- true/false

params.attSchedule
- 意味：スケジュール本体（struct配列）
- フィールド
  - tStart：開始時刻（datetime, UTC）
  - tEnd：終了時刻（datetime, UTC）
  - mode：fixed/keyframes/nadir/sun/velocity のいずれか
  - opt：その区間で使う attOpt（または keyframes 用の差し替え設定）

区間外の扱い
- 実装方針例
  - 最後の区間を採用する
  - 先頭の区間を採用する
  - 既定の params.attMode/params.attOpt にフォールバックする
- 実装に合わせて運用を決める

---

## 11. STLモデル（衛星形状）表示関連

params.useSatModel
- 意味：STLモデル表示を使うか
- true/false

params.satModelFile
- 意味：STLファイル名またはパス
- 例：satellite.stl
- 注意：カレントフォルダに置くか、フルパス指定する

params.satModelUnit
- 意味：STLモデルの単位（STL自体は単位情報を持たないため手動指定）
- 値：mm / m / km
- 例：mm の場合は 1e-6 で km に変換

params.satModelDesiredSize_km_main
- 意味：右図メインの衛星モデル表示サイズ（目標）
- 単位：km
- 注意：大きくしすぎると地球や軌道線を覆う。視認性とのバランス調整が必要

params.satModelDesiredSize_km_inset
- 意味：インセットの衛星モデル表示サイズ（目標）
- 単位：km
- 用途：インセットは確認用なので大きめで問題ない

params.fallbackMarkerColor
- 意味：STLが読めない場合に使う点表示の色（RGB）

---

## 12. 推奨設定例

インセットを小さめ＋右上詰め（被り軽減）
```matlab
params.insetSizeFactor = 0.15;
params.insetMargin = 0.002;
````

時間刻みを1分に

```matlab
params.timeStepMin = 1;
```

Sun pointing（+Xをsun、+Zをnadir、X軸まわりに30度ロール）

```matlab
params.attMode = "sun";
params.attOpt = struct('xAxis',"sun",'zAxis',"nadir",'rollDeg',30,'rollAxis',"x");
```

区間スケジュール（例：0-6h nadir、6-12h sun、12-24h velocity）

```matlab
params.useAttSchedule = true;

params.attSchedule(1) = struct( ...
 'tStart', datetime(params.baseYear,1,1,0,0,0,'TimeZone','UTC'), ...
 'tEnd', datetime(params.baseYear,1,1,6,0,0,'TimeZone','UTC'), ...
 'mode', "nadir", ...
 'opt', struct('xAxis',"vel",'zAxis',"nadir",'rollDeg',15,'rollAxis',"z") );

params.attSchedule(2) = struct( ...
 'tStart', datetime(params.baseYear,1,1,6,0,0,'TimeZone','UTC'), ...
 'tEnd', datetime(params.baseYear,1,1,12,0,0,'TimeZone','UTC'), ...
 'mode', "sun", ...
 'opt', struct('xAxis',"sun",'zAxis',"nadir",'rollDeg',30,'rollAxis',"x") );

params.attSchedule(3) = struct( ...
 'tStart', datetime(params.baseYear,1,1,12,0,0,'TimeZone','UTC'), ...
 'tEnd', datetime(params.baseYear,1,2,0,0,0,'TimeZone','UTC'), ...
 'mode', "velocity", ...
 'opt', struct('xAxis',"vel",'zAxis',"nadir",'rollDeg',0,'rollAxis',"x") );
```

***

## 13. トラブルシューティング（よくある症状）

STLが点表示になる

*   STL読み込みに失敗している可能性がある
*   satellite.stl の配置（カレントフォルダ）やパス指定を確認する
*   satModelUnit が合っているか確認する

トライアドが見えない

*   showInsetTriad が true か確認する
*   triadAutoFrac と insetHalfWidth\_km の組合せで矢印が極端に短くなっていないか確認する
*   インセット作成後に quiver を作っているか確認する（作成順序が重要）

再生（Play）が動かない

*   timer が生成されているか確認する
*   CloseRequestFcn で timer を停止・削除しているか確認する
*   isAutoUpdating の扱いが誤っていると更新が止まる場合がある
