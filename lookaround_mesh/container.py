"""MCP4 容器 → LZFSE → CHNK 小区画、と各区画の見出し。"""
import struct

import numpy as np

# MCP4 の記録の種類
MESH, FACE_INDEX, FACE_IMAGE, INDEX_TABLE = 1, 2, 3, 0
# メッシュ内の CHNK の種類
CHNK_TRANSFORM, CHNK_CLERS, CHNK_FACE_ID, CHNK_VERTICES, CHNK_HEADER_COPY, CHNK_PARAMS = 0, 4, 5, 9, 10, 11
CHNK_HEADER = 22  # 'CHNK', u16 版, u16 種類, ..., u32 本体の大きさ


def mcp4_records(b):
    """[(種類, 位置, 大きさ)]。0x12 から 10 バイトずつ、ファイルをすき間なく覆う。"""
    if b[:4] != b"MCP4":
        raise ValueError(f"not an MCP4 container: {b[:4]!r}")
    recs, o, expect = [], 0x12, None
    while o + 10 <= len(b):
        t, off, size = struct.unpack_from("<HII", b, o)
        if (expect is not None and off != expect) or off + size > len(b) or t > 64:
            break
        recs.append((t, off, size))
        expect, o = off + size, o + 10
    return recs


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
