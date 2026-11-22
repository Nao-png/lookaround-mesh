"""MCP4 容器 → LZFSE → CHNK 小区画、と各区画の見出し。"""
import struct

import numpy as np

# MCP4 の記録の種類
MESH, FACE_INDEX, FACE_IMAGE, INDEX_TABLE, PANO_INFO = 1, 2, 3, 0, 5
# メッシュ内の CHNK の種類
CHNK_TRANSFORM, CHNK_CLERS, CHNK_FACE_ID, CHNK_VERTICES, CHNK_PROJECTION, CHNK_PARAMS = 0, 4, 5, 9, 10, 11
CHNK_HEADER_COPY = CHNK_PROJECTION  # 旧名
CHNK_HEADER = 22  # 'CHNK', u16 版, u16 種類, ..., u32 本体の大きさ
TABLE = 0x08      # 記録の表の始まり。記録の数は 0x06 の u16


def mcp4_records(b):
    """[(種類, 位置, 大きさ)]。0x08 から 10 バイトずつ、0x06 の u16 の数だけ。
    記録はファイルの残り（表の直後から）をすき間なく、順に覆う。"""
    if b[:4] != b"MCP4":
        raise ValueError(f"not an MCP4 container: {b[:4]!r}")
    n = struct.unpack_from("<H", b, 6)[0]
    recs = [struct.unpack_from("<HII", b, TABLE + 10 * i) for i in range(n)]
    expect = TABLE + 10 * n
    for t, off, size in recs:
        if off != expect or off + size > len(b):
            raise ValueError(f"MCP4 record ({t}, {off:#x}, {size}) does not follow at {expect:#x}")
        expect = off + size
    return recs


def pano_info(b):
    """MCP4 の種類 5（224 バイト）。パノラマ id、ECEF 位置、方位ごとの表。

    u16 0, u64 id, f64×3 ECEF, f32 0.01, u16 180, u32 127, 180 バイト（値 = バイト − 127）。
    180 個の値は、先頭が最初の値で残りが差。累積和 profile[i]（0〜100）が、面 0 の中心から
    右回りに方位 2i° の値（意味は未確定、docs/FORMAT.md）。
    """
    (_, off, size), = [r for r in mcp4_records(b) if r[0] == PANO_INFO]
    r = b[off:off + size]
    panoid, = struct.unpack_from("<Q", r, 2)
    ecef = np.array(struct.unpack_from("<3d", r, 10))
    scale, n, bias = struct.unpack_from("<fHI", r, 34)
    values = np.frombuffer(r, np.uint8, n, 44).astype(np.int32) - bias
    return {"panoid": panoid, "ecef": ecef, "scale": scale, "values": values, "profile": np.cumsum(values)}


def meshes(b):
    """各面のメッシュ（展開済み）。種類 1 = u8 3、u32 展開後の大きさ、LZFSE。"""
    import liblzfse

    out = []
    for t, off, size in mcp4_records(b):
        if t == MESH:
            d = b[off:off + size]
            raw = liblzfse.decompress(d[5:])
            want = struct.unpack_from("<I", d, 1)[0]
            if len(raw) != want:
                raise ValueError(f"LZFSE gave {len(raw)} bytes, header says {want}")
            out.append(raw)
    return out


def face_images(b):
    """各面の低解像度画像（HEIC のバイト列）。記録の先頭 1 バイトは HEIC の前置き。"""
    return [b[off + 1:off + size] for t, off, size in mcp4_records(b) if t == FACE_IMAGE]


def chnks(m):
    """{種類: 本体}。先頭の種類 0 だけは見出しごと返す（行列の位置が見出し基準のため）。"""
    offs = [i for i in range(len(m)) if m[i:i + 4] == b"CHNK"] + [len(m)]
    out = {}
    for a, b in zip(offs, offs[1:]):
        t = struct.unpack_from("<H", m, a + 6)[0]
        out[t] = m[a:b] if t == CHNK_TRANSFORM else m[a + CHNK_HEADER:b]
    return out


def header(m):
    """(局所→地心座標 ECEF の 4x4 行列, 頂点数 V, 三角形数 F)。"""
    mat = np.array(struct.unpack_from("<16d", m, 0x1e)).reshape(4, 4)
    nv, nf = struct.unpack_from("<II", m, 0x9e)
    return mat, nv, nf


def projection(m):
    """面の投影 (su, sv, ou, ov)。種類 10 の末尾の f64×4。

    角度 [rad] → 面の画像の正規化座標 [-1, 1]: (u·su + ou, v·sv + ov)。被覆タイルの
    レンズ投影と同じ値で、su = 2/fov_s、sv = 2/fov_h、ou = cx·su、ov = cy·sv。
    """
    return struct.unpack_from("<4d", chnks(m)[CHNK_PROJECTION], 142)
