"""Apple Look Around の mt/7 ファイル（パノラマごとのメッシュ）の解読。

MCP4 容器 → LZFSE → CHNK 小区画 → 位置の行列、Edgebreaker の接続情報、
頂点の値と、その座標への組み立て（docs/FORMAT.md）。
"""
from .container import chnks, face_images, header, mcp4_records, meshes, pano_info, projection
from .edgebreaker import clers, decode, first_touch, split_lengths
from .vertices import apple_ids, decode_vertices, to_ecef, to_local, vertex_values


def face_mesh(m):
    """1 面分: (三角形の頂点番号の組, (V, 3) の整数座標 (u, v, w))。"""
    s, b = clers(m)
    tris, _, _ = decode(s, b)
    return tris, decode_vertices(m, tris, b)


__all__ = ["mcp4_records", "meshes", "chnks", "header", "face_images", "pano_info", "projection",
           "clers", "split_lengths", "decode", "first_touch",
           "vertex_values", "apple_ids", "decode_vertices", "to_local", "to_ecef",
           "face_mesh"]
