# experiments/

頂点の座標を組み立てる規則を探すのに使ったスクリプト。結果と、物差しの落とし穴は [../docs/HYPOTHESES.md](../docs/HYPOTHESES.md)。

- `harness.py` — サンプルのメッシュをまとめて読み込み、「発散の位置」で候補を採点する。**冒頭から発散する候補を見逃す**（HYPOTHESES.md の 4）。
- `probe.py` — 1 か所だけ予測の規則を差し替えて組み立て直し、「浮いた頂点」が最初に現れる位置を測る。**領域ごとまとめてずれると見逃す**（HYPOTHESES.md の 6）。

どちらも `samples/` を読む。
