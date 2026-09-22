"""Apple Look Around の mt/7 ファイル（パノラマごとのメッシュ）の解読。

解けている層: MCP4 容器、LZFSE、CHNK 小区画、位置の行列、残差の詰め方、
Edgebreaker の接続情報。解けていない層: 頂点の予測規則と座標の意味
（docs/FORMAT.md、docs/HYPOTHESES.md）。
"""
from .container import chnks, face_images, header, mcp4_records, meshes
from .edgebreaker import clers, decode, first_touch, split_lengths
from .vertices import vertex_values

__all__ = ["mcp4_records", "meshes", "chnks", "header", "face_images",
           "clers", "split_lengths", "decode", "first_touch", "vertex_values"]
