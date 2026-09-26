# samples/ について

渋谷スクランブル交差点の Look Around パノラマ 8 枚（`neighbours.json` に id・撮影日時・中心からの距離などの一覧）。

mt/7 ファイル（`pano*_mt7.bin`）は Apple Look Around のサーバーから取得するデータで、**Apple の著作物**です。そのためリポジトリには入れていません。`python scripts/fetch_samples.py` で取ってください（`.gitignore` 済み）。

`docs/images/` の画像（README の写真）は、解読の結果を示すために Look Around のデータ（面の画像とメッシュ）から作ったもので、Apple の画像を含みます。このリポジトリのライセンス（MIT）の対象ではありません。

Apple Maps のデータをこの形で取り出すことは、Apple の利用規約に抵触しうる点に注意してください。取り出したデータは再配布しないでください。

---

The 8 Look Around panoramas at the Shibuya scramble crossing (listed with ids, capture times and distances in `neighbours.json`).

The mt/7 files (`pano*_mt7.bin`) are Apple's copyrighted data fetched from Apple's Look Around servers, so they are not in the repository. Fetch them with `python scripts/fetch_samples.py` (they are git-ignored).

The pictures in `docs/images/` were made from Look Around data (face images and meshes) to illustrate the results. They contain Apple imagery and are not covered by this repository's licence (MIT).

Extracting Apple Maps data this way may violate Apple's terms of service. Do not redistribute the data you extract.
