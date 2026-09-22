"""mt/7 の取得。認証は streetlevel の Authenticator をそのまま使う。"""
import requests
from streetlevel.lookaround.auth import Authenticator

ENDPOINT = "https://gspe72-ssl.ls.apple.com/mnn_us/"
_auth = None


def mt_url(panoid, build_id):
    global _auth
    _auth = _auth or Authenticator()
    s = str(panoid).zfill(20)
    path = "/".join(s[i:i + 4] for i in range(0, 20, 4))
    return _auth.authenticate_url(f"{ENDPOINT}{path}/{str(build_id).zfill(10)}/mt/7")


def fetch_mt(panoid, build_id, timeout=60):
    r = requests.get(mt_url(panoid, build_id), timeout=timeout)
    r.raise_for_status()
    return r.content
