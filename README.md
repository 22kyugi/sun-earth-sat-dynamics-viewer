# sun-earth-sat-dynamics-viewer

太陽中心（ヘリオセントリック）視点で **地球の公転**と **衛星の地球周回**を可視化する MATLAB スクリプトです。  
UI（スライダー）で **月 / 日 / 時刻**を直感的に操作し、任意時刻の **太陽方向・昼夜境界（ターミネータ）・衛星位置・軌道面**を表示できます。

> English (short): A MATLAB visualization tool for Sun-centered Earth orbit and Earth-centered satellite motion with UI sliders (month/day/hour). It renders satellite CAD models from STL and shows the instantaneous orbit plane.

---

## Features

- **Sun-centered view (left)**
  - 地球の公転軌道（1年分）
  - 地球の移動（テクスチャ付き地球、昼夜陰影、ターミネータ）
  - 春分・夏至・秋分・冬至のマーカー表示（近似日付）
  - 太陽の橙色マーカー表示（原点）

- **Earth-centered view (right)**
  - 衛星の現在位置（STLモデル表示 / フォールバック点）
  - 現在時刻の **瞬時軌道（軌道線）** のみを表示
  - **軌道面（半透明ディスク）** + 法線ベクトル矢印
  - 簡易 **地球影（食）判定**（円柱影近似）

- **UI**
  - メインスライダー：月（`2000+n年1月`〜`2000+n+1年1月`）
  - サブスライダー：日にち（1〜月末）
  - サブスライダー：時刻（0〜24h）

---

## Demo / Screenshot

> ここにスクリーンショットを追加してください  
例：`docs/screenshot.png`

```text
docs/
  screenshot.png
