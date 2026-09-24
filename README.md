# lookaround-mesh

Apple Look Around のパノラマごとに配信されている 3D メッシュ（`mt/7` ファイル）の解読記録とコード。

Look Around のパノラマには、6 つの面（経度で 4 つに切った帯と、上下の帽子）それぞれについて小さなメッシュが付いている。形式は非公開で、Look Around を解析してきた [sk-zk/lookaround](https://github.com/sk-zk/lookaround) でも「容器を区画に切り分けるところまで、区画の中身のパーサはまだ無い」とされていた。ここでは頂点の座標まで、全層を解いた。

## 状況

| 層 | 状態 | 裏付け |
|---|---|---|
| 容器 `MCP4`（記録の表） | ✅ 解読 | 8/8 のファイルで記録がすき間なく全体を覆う |
| 圧縮（LZFSE） | ✅ 解読 | 48/48 のメッシュが宣言どおりの大きさに展開 |
| 小区画 `CHNK` の並び | ✅ 解読 | 種類 0・4・5・9・10・11 |
| 位置の行列（局所 → 地心座標 ECEF） | ✅ 解読 | 回転が正規直交、パノラマ間の距離が GPS と 0.2 m 以内（8/8） |
| 接続情報（Edgebreaker の CLERS 文字列） | ✅ 解読 | 48/48 で三角形数・頂点数一致、向きの一貫した円盤状メッシュに復元 |
| 頂点の残差の詰め方（1 バイト＋例外表） | ✅ 解読 | 48/48 で例外表の長さが 0xff の数と一致 |
| 頂点の並び（軸ごと、Apple の頂点番号順） | ✅ 解読 | 最初の三角形の 3 頂点（生の値）の位置を 48/48 で当てる |
| 座標の組み立て（面の走査＋平行四辺形予測） | ✅ 解読 | 48/48 で全頂点が面の範囲に収まり、範囲が 8 枚で共通 |
| 座標の意味（経度・緯度・距離の逆数） | ✅ 解読 | 地面が平面でカメラ高 2.30〜2.32 m（8/8）。翌日の別走行のメッシュと地心座標で重なる（中央値 0.25 m、対照 4.5 m） |

つまり、各メッシュの三角形と頂点の地心座標（ECEF）がすべて得られる。

形式の詳細は [docs/FORMAT.md](docs/FORMAT.md)。

## 使い方

```bash
pip install -r requirements.txt
python scripts/check.py            # 解けている層の主張をすべて samples/ で検証
python -m pytest tests -q
python scripts/fetch_samples.py    # samples/ を取り直す（Apple のサーバーから）
```

`scripts/check.py` の出力例（すべて OK になる）:

```
OK MCP4: 記録がファイルをすき間なく覆う                       8/8
OK 最初のループの長さ（文字列から計算）= V - Cの数          48/48
OK 復元: どの辺も 2 枚以下、向きが一貫                      48/48
OK 頂点: 例外表の長さ = 本体の 0xff の数                    48/48
OK 位置: パノラマ間の距離が GPS と 0.2 m 以内                8/8
...
```

コードから使う:

```python
import lookaround_mesh as lm

b = open("samples/pano0_mt7.bin", "rb").read()
for raw in lm.meshes(b):                 # 立方体の 6 面
    mat, nv, nf = lm.header(raw)         # 局所→ECEF の 4x4 行列、頂点数、三角形数
    s, boundary = lm.clers(raw)          # CLERS 文字列と境界の頂点数
    tris, gates, _ = lm.decode(s, boundary)   # 三角形（頂点番号の組）
    P = lm.decode_vertices(raw, tris, boundary)  # (V, 3): 経度・緯度 [mrad]、距離の逆数 w
    xyz = lm.to_ecef(raw, P)             # 地心座標 [m]（空の頂点 w = 0 は inf）
    # まとめて: tris, P = lm.face_mesh(raw)
```

OBJ に書き出す（パノラマの位置を原点にした東・北・上）:

```bash
python scripts/export_obj.py samples/pano0_mt7.bin pano0.obj --max-edge 3
```

色付きの点群（PLY）に書き出す。三角形の面に点をまき、面の画像から色を引き、複数のパノラマを 1 つにまとめる:

```bash
python scripts/export_ply.py samples/pano*_mt7.bin cloud.ply --spacing 0.1 --max-range 40
```

ほかの場所のパノラマを取ってきて、点群にし、ブラウザで見る:

```bash
python scripts/fetch_area.py 36.5722 136.6680 --radius 250 --out out/kanazawa
python scripts/export_ply.py out/kanazawa out/kanazawa.ply --spacing 0.1 --max-range 15 --jobs 12
python scripts/export_web.py out/kanazawa.ply out/web
python scripts/serve_web.py out/web 8000        # ブラウザで http://localhost:8000
```

メッシュのまま、面の画像を貼ってブラウザで見る:

```bash
python scripts/export_web_mesh.py out/kanazawa out/web_mesh --jobs 12
python scripts/serve_web.py out/web_mesh 8000
```

場所ごとに、それをいちばん正確に写すパノラマ（面上の 1 画素の大きさと奥行きの刻みで比べる）のメッシュだけを使う。

渋谷のように撮影の多い場所では、年をそろえて間引いてから取る（建て替えが多く、年をまたぐと建物が食い違う）。ビルが高いので距離の上限を広げる:

```bash
python scripts/fetch_area.py 35.6595 139.7005 --radius 250 --date 2026 --min-spacing 3 --out out/shibuya
python scripts/export_web_mesh.py out/shibuya out/web_shibuya --jobs 10 --max-range 50 --min-facing 0.15 --sky-range 30 --max-edge 4
```

奥行き（距離の逆数 w）は整数なので、距離 r での刻みは約 r² / 1000 m（50 m 先で 2.5 m）。高いビルの上部は粗く、破片のようになる。

## 構成

```
lookaround_mesh/    解けている層のライブラリ
  container.py      MCP4、LZFSE、CHNK、行列と頂点数
  edgebreaker.py    CLERS 文字列から三角形を復元
  vertices.py       頂点の値、頂点番号、座標の組み立て、ECEF への変換
  fetch.py          mt/7 の取得
scripts/check.py    すべての主張の検証
scripts/fetch_samples.py
scripts/export_obj.py    OBJ に書き出す
scripts/export_ply.py    色付きの点群を PLY に書き出す
scripts/fetch_area.py    ある地点のまわりの mt/7 をまとめて取る
scripts/export_web.py    点群をブラウザ用（three.js）に変換する
scripts/export_web_mesh.py  メッシュを画像付きでブラウザ用に変換する
scripts/serve_web.py     その出力を手元で配る
experiments/        座標復元の探索で使った初期のスクリプト（docs/HYPOTHESES.md）
samples/            渋谷スクランブル交差点のパノラマ 8 枚の mt/7（下記の注意）
docs/
```

## サンプル

渋谷スクランブル交差点（35.6595, 139.7005）の Look Around パノラマ 8 枚:

- 中心（id 14755389547370281509、2026-04-21 の車載撮影）と、同じ走行で 4.5 m 間隔の 5 枚
- 翌日（2026-04-22）の別走行で、中心から 1.7 m と 3.3 m の 2 枚

別走行の 2 枚は同じ建物を別の位置から見ているので、座標の解読が正しければ、両者のメッシュは地心座標で重なるはず。これを座標復元の判定基準にし、満たした（docs/FORMAT.md の「検証」）。

## 注意: データの扱い

`samples/` の `.bin` は Apple の著作物で、このレポジトリのライセンスの対象外（[samples/NOTICE.md](samples/NOTICE.md)）。このレポジトリは非公開を前提にしている。**公開に切り替えるなら、先に `samples/` の `.bin` を履歴ごと取り除くこと。** Apple Maps のデータをこの形で取り出すことは Apple の利用規約に抵触しうる。

## 先行研究

- [sk-zk/streetlevel](https://github.com/sk-zk/streetlevel) — Look Around の取得と認証。`fetch.py` はその `Authenticator` を使う。
- [sk-zk/lookaround](https://github.com/sk-zk/lookaround) — `mt/7` の取得と MCP4 の区画分け。
- [retroplasma/flyover-reverse-engineering](https://github.com/retroplasma/flyover-reverse-engineering) — Apple Flyover の C3M 形式の Edgebreaker 変種。頂点番号の振り方と面の走査は同じ系統で、座標の組み立てを解く手がかりになった。

## ライセンス

コードは MIT（[LICENSE](LICENSE)）。`samples/` のデータは含まない。
