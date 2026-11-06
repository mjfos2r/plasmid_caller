"""Tests for Cloudflare Access middleware."""

import base64
import json
import time
from types import SimpleNamespace

import pytest

jwt = pytest.importorskip("jwt")
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import rsa
from fastapi import FastAPI, Request
from fastapi.testclient import TestClient
from jwt import PyJWKClient
from jwt.algorithms import RSAAlgorithm
from jwt.exceptions import PyJWKClientConnectionError

from plasmid_caller.web.auth import AccessConfig, CloudflareAccessMiddleware, install_access_auth

CFG = AccessConfig(team_domain="testteam.cloudflareaccess.com", audience="test-aud")

GOOD_KEY = rsa.generate_private_key(public_exponent=65537, key_size=2048)
EVIL_KEY = rsa.generate_private_key(public_exponent=65537, key_size=2048)


class FakeJWKS:
    """Always returns the legitimate public key, as the real JWKS would."""

    def get_signing_key_from_jwt(self, token):
        return SimpleNamespace(key=GOOD_KEY.public_key())


def make_token(key=GOOD_KEY, kid="k1", **overrides):
    now = int(time.time())
    claims = {
        "aud": [CFG.audience],
        "iss": CFG.issuer,
        "iat": now,
        "exp": now + 600,
        "email": "pi@example.org",
    }
    claims.update(overrides)
    claims = {k: v for k, v in claims.items() if v is not None}
    return jwt.encode(claims, key, algorithm="RS256", headers={"kid": kid})


def client_with(jwks):
    app = FastAPI()

    @app.get("/whoami")
    def whoami(request: Request):
        return {"email": request.state.access_claims.get("email")}

    app.add_middleware(CloudflareAccessMiddleware, config=CFG, jwks_client=jwks)
    return TestClient(app)


@pytest.fixture
def client():
    return client_with(FakeJWKS())


def real_jwks(*keys):
    """A real PyJWKClient, with only the network fetch replaced."""
    jwks = PyJWKClient(CFG.certs_url, cache_keys=True, lifespan=3600)
    entries = []
    for kid, key in keys:
        entry = json.loads(RSAAlgorithm.to_jwk(key.public_key()))
        entries.append({**entry, "kid": kid, "use": "sig", "alg": "RS256"})
    jwks.fetch_data = lambda: {"keys": entries}
    return jwks


def call(client, token=None):
    headers = {"Cf-Access-Jwt-Assertion": token} if token else {}
    return client.get("/whoami", headers=headers)


def test_valid_token_passes_and_exposes_email(client):
    r = call(client, make_token())
    assert r.status_code == 200
    assert r.json() == {"email": "pi@example.org"}


def test_missing_header_rejected(client):
    assert call(client).status_code == 403


def test_garbage_token_rejected(client):
    assert call(client, "not.a.jwt").status_code == 403


def test_forged_signature_rejected(client):
    assert call(client, make_token(key=EVIL_KEY)).status_code == 403


def test_wrong_audience_rejected(client):
    assert call(client, make_token(aud=["some-other-app"])).status_code == 403


def test_wrong_issuer_rejected(client):
    assert call(client, make_token(iss="https://evil.cloudflareaccess.com")).status_code == 403


def test_expired_rejected(client):
    past = int(time.time()) - 3600
    assert call(client, make_token(iat=past - 600, exp=past)).status_code == 403


def test_missing_exp_rejected(client):
    assert call(client, make_token(exp=None)).status_code == 403


def test_alg_none_rejected(client):
    b64 = lambda d: base64.urlsafe_b64encode(json.dumps(d).encode()).rstrip(b"=").decode()
    now = int(time.time())
    payload = {"aud": [CFG.audience], "iss": CFG.issuer, "iat": now, "exp": now + 600}
    token = f"{b64({'alg': 'none', 'typ': 'JWT'})}.{b64(payload)}."
    assert call(client, token).status_code == 403


def test_refuses_to_start_without_config(monkeypatch):
    monkeypatch.delenv("CF_TEAM_DOMAIN", raising=False)
    monkeypatch.delenv("CF_ACCESS_AUD", raising=False)
    monkeypatch.delenv("ALLOW_INSECURE_DEV", raising=False)
    with pytest.raises(RuntimeError):
        install_access_auth(FastAPI())


def test_insecure_dev_blocked_on_cloud_run(monkeypatch):
    monkeypatch.setenv("ALLOW_INSECURE_DEV", "1")
    monkeypatch.setenv("K_SERVICE", "plasmid-portal")
    with pytest.raises(RuntimeError):
        install_access_auth(FastAPI())


def test_real_client_accepts_the_published_key():
    client = client_with(real_jwks(("k1", GOOD_KEY)))
    assert call(client, make_token()).status_code == 200


def test_real_client_rejects_a_token_signed_by_another_key_under_a_known_kid():
    client = client_with(real_jwks(("k1", GOOD_KEY)))
    assert call(client, make_token(key=EVIL_KEY)).status_code == 403


def test_real_client_rejects_an_unknown_kid():
    client = client_with(real_jwks(("k1", GOOD_KEY)))
    assert call(client, make_token(kid="unknown")).status_code == 403


def test_real_client_picks_up_a_rotated_key():
    jwks = real_jwks(("k1", GOOD_KEY))
    client = client_with(jwks)
    assert call(client, make_token()).status_code == 200
    rotated = real_jwks(("k1", GOOD_KEY), ("k2", EVIL_KEY))
    jwks.fetch_data = rotated.fetch_data
    assert call(client, make_token(key=EVIL_KEY, kid="k2")).status_code == 200


def test_public_key_used_as_hmac_secret_rejected():
    pem = GOOD_KEY.public_key().public_bytes(serialization.Encoding.PEM,
                                             serialization.PublicFormat.SubjectPublicKeyInfo)
    b64 = lambda d: base64.urlsafe_b64encode(d).rstrip(b"=").decode()
    now = int(time.time())
    head = b64(json.dumps({"alg": "HS256", "kid": "k1"}).encode())
    body = b64(json.dumps({"aud": [CFG.audience], "iss": CFG.issuer, "iat": now, "exp": now + 600}).encode())
    import hashlib, hmac
    signature = b64(hmac.new(pem, f"{head}.{body}".encode(), hashlib.sha256).digest())
    client = client_with(real_jwks(("k1", GOOD_KEY)))
    assert call(client, f"{head}.{body}.{signature}").status_code == 403


@pytest.mark.parametrize("token", [".", "..", "a.b", "a.b.c.d", "e30.e30.", "eyJhbGciOiJSUzI1NiJ9..", "x" * 5000])
def test_malformed_tokens_are_403_not_500(token):
    client = TestClient(client_with(real_jwks(("k1", GOOD_KEY))).app, raise_server_exceptions=False)
    assert call(client, token).status_code == 403


def test_jwks_outage_is_503():
    class Down:
        def get_signing_key_from_jwt(self, token):
            raise PyJWKClientConnectionError("down")

    assert call(client_with(Down()), make_token()).status_code == 503


def test_exempt_prefix_bypasses_auth_but_only_that_prefix():
    app = FastAPI()

    @app.get("/api/public/anything")
    def public_anything():
        return {"ok": True}

    @app.get("/api/public-not-really/x")
    def not_really():
        return {"ok": True}

    app.add_middleware(CloudflareAccessMiddleware, config=CFG, jwks_client=FakeJWKS(),
                       exempt_prefixes=frozenset({"/api/public/"}))
    client = TestClient(app)
    assert client.get("/api/public/anything").status_code == 200
    assert client.get("/api/public-not-really/x").status_code == 403


def test_insecure_dev_allowed_locally(monkeypatch):
    monkeypatch.setenv("ALLOW_INSECURE_DEV", "1")
    monkeypatch.delenv("K_SERVICE", raising=False)
    app = FastAPI()
    install_access_auth(app)
    assert not app.user_middleware
