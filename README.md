# lookaround-mesh

Apple Look Around のパノラマごとに配信されている 3D メッシュ（`mt/7` ファイル）の解読記録とコード。

Look Around のパノラマには、立方体の 6 面それぞれについて小さなメッシュが付いている。形式は非公開で、Look Around を解析してきた [sk-zk/lookaround](https://github.com/sk-zk/lookaround) でも「容器を区画に切り分けるところまで、区画の中身のパーサはまだ無い」とされていた。ここではその先を解いた。

## 状況

| 層 | 状態 | 裏付け |
|---|---|---|
| 容器 `MCP4`（記録の表） | ✅ 解読 | 8/8 のファイルで記録がすき間なく全体を覆う |
| 圧縮（LZFSE） | ✅ 解読 | 48/48 のメッシュが宣言どおりの大きさに展開 |
| 小区画 `CHNK` の並び | ✅ 解読 | 種類 0・4・5・9・10・11 |
| 位置の行列（局所 → 地心座標 ECEF） | ✅ 解読 | 回転が正規直交、パノラマ間の距離が GPS と 0.2 m 以内（8/8） |
| 接続情報（Edgebreaker の CLERS 文字列） | ✅ 解読 | 48/48 で三角形数・頂点数一致、向きの一貫した円盤状メッシュに復元 |
| 頂点の残差の詰め方（1 バイト＋例外表） | ✅ 解読 | 48/48 で例外表の長さが 0xff の数と一致 |
| 頂点の並び順 | △ 一部 | 境界は「走査で初めて触れた順」でまとまるが、全体は未確定 |
| **予測の規則と座標の意味** | ❌ 未解明 | [docs/HYPOTHESES.md](docs/HYPOTHESES.md) に試して外れた仮説を記録 |

つまり「どの頂点とどの頂点が三角形を作るか」と「各メッシュが地球上のどこにあるか」は分かる。「各頂点がどこにあるか」はまだ分からない。

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
    values = lm.vertex_values(raw)       # 頂点の残差 V×3 個（座標への組み立ては未解明）
```

## 構成

```
lookaround_mesh/    解けている層のライブラリ
  container.py      MCP4、LZFSE、CHNK、行列と頂点数
  edgebreaker.py    CLERS 文字列から三角形を復元
  vertices.py       頂点の残差の取り出し（＋未完成の座標復元の試み）
  fetch.py          mt/7 の取得
scripts/check.py    すべての主張の検証
scripts/fetch_samples.py
experiments/        座標復元の探索に使ったスクリプト（結果は docs/HYPOTHESES.md）
samples/            渋谷スクランブル交差点のパノラマ 8 枚の mt/7（下記の注意）
docs/
```

## サンプル

渋谷スクランブル交差点（35.6595, 139.7005）の Look Around パノラマ 8 枚:

- 中心（id 14755389547370281509、2026-04-21 の車載撮影）と、同じ走行で 4.5 m 間隔の 5 枚
- 翌日（2026-04-22）の別走行で、中心から 1.7 m と 3.3 m の 2 枚

別走行の 2 枚は同じ建物を別の位置から見ているので、座標の解読が正しければ、両者のメッシュは地心座標で重なるはず。これが座標復元の判定基準になる。

## 注意: データの扱い

`samples/` の `.bin` は Apple の著作物で、このレポジトリのライセンスの対象外（[samples/NOTICE.md](samples/NOTICE.md)）。このレポジトリは非公開を前提にしている。**公開に切り替えるなら、先に `samples/` の `.bin` を履歴ごと取り除くこと。** Apple Maps のデータをこの形で取り出すことは Apple の利用規約に抵触しうる。

## 先行研究

- [sk-zk/streetlevel](https://github.com/sk-zk/streetlevel) — Look Around の取得と認証。`fetch.py` はその `Authenticator` を使う。
- [sk-zk/lookaround](https://github.com/sk-zk/lookaround) — `mt/7` の取得と MCP4 の区画分け。

## ライセンス

コードは MIT（[LICENSE](LICENSE)）。`samples/` のデータは含まない。
