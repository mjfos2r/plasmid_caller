"""Cloudflare Access JWT verification for the plasmid portal.

Pure ASGI middleware. Every HTTP/WebSocket request must carry a valid
`Cf-Access-Jwt-Assertion` header issued by your Cloudflare Access team for
this application (checked: RS256 signature against the team JWKS, aud, iss,
exp, iat). Anything else is rejected before it reaches a route.

Fails closed: the app refuses to start if CF_TEAM_DOMAIN / CF_ACCESS_AUD are
unset. Local development can opt out with ALLOW_INSECURE_DEV=1, which is
itself refused when running on Cloud Run (K_SERVICE is set there).

Env:
    CF_TEAM_DOMAIN   e.g. yourteam.cloudflareaccess.com
    CF_ACCESS_AUD    Application Audience (AUD) tag from the Access app

Usage:
    from plasmid_caller.web.auth import install_access_auth
    app = FastAPI()
    install_access_auth(app)
    # in routes: request.state.access_claims["email"]
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass
from typing import Any, Protocol

import anyio
import jwt
from jwt import PyJWKClient
from jwt.exceptions import PyJWKClientConnectionError, PyJWTError
from starlette.types import ASGIApp, Receive, Scope, Send

log = logging.getLogger("portal.auth")

_HEADER = b"cf-access-jwt-assertion"
_REQUIRED_CLAIMS = ["exp", "iat", "iss", "aud"]


@dataclass(frozen=True)
class AccessConfig:
    team_domain: str
    audience: str

    @property
    def issuer(self) -> str:
        return f"https://{self.team_domain}"

    @property
    def certs_url(self) -> str:
        return f"{self.issuer}/cdn-cgi/access/certs"

    @classmethod
    def from_env(cls) -> "AccessConfig":
        team = (
            os.environ.get("CF_TEAM_DOMAIN", "")
            .strip()
            .removeprefix("https://")
            .rstrip("/")
        )
        aud = os.environ.get("CF_ACCESS_AUD", "").strip()
        if not team or not aud:
            raise RuntimeError(
                "CF_TEAM_DOMAIN and CF_ACCESS_AUD must be set. "
                "Refusing to start without authentication."
            )
        return cls(team_domain=team, audience=aud)


class SigningKeySource(Protocol):
    def get_signing_key_from_jwt(self, token: str) -> Any: ...


class CloudflareAccessMiddleware:
    def __init__(
        self,
        app: ASGIApp,
        config: AccessConfig,
        jwks_client: SigningKeySource | None = None,
        exempt_paths: frozenset[str] = frozenset(),
        exempt_prefixes: frozenset[str] = frozenset(),
        leeway_seconds: int = 30,
    ) -> None:
        self.app = app
        self.config = config
        # PyJWKClient caches keys and refetches on unknown `kid` (key rotation).
        self.jwks = jwks_client or PyJWKClient(
            config.certs_url, cache_keys=True, lifespan=3600
        )
        self.exempt_paths = exempt_paths
        self.exempt_prefixes = exempt_prefixes
        self.leeway = leeway_seconds

    def _verify(self, token: str) -> dict[str, Any]:
        signing_key = self.jwks.get_signing_key_from_jwt(token)
        return jwt.decode(
            token,
            signing_key.key,
            algorithms=["RS256"],  # pinned; never trust the token's own alg
            audience=self.config.audience,
            issuer=self.config.issuer,
            leeway=self.leeway,
            options={"require": _REQUIRED_CLAIMS},
        )

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        if scope["type"] not in ("http", "websocket"):
            await self.app(scope, receive, send)  # lifespan etc.
            return

        path = scope.get("path")
        if path in self.exempt_paths or any(path.startswith(prefix) for prefix in self.exempt_prefixes):
            await self.app(scope, receive, send)
            return

        token = None
        for name, value in scope.get("headers", []):
            if name == _HEADER:
                token = value.decode("latin-1")
                break

        if not token:
            await _deny(scope, send, 403, "missing access token")
            return

        try:
            # JWKS fetch is blocking I/O; keep it off the event loop.
            claims = await anyio.to_thread.run_sync(self._verify, token)
        except PyJWKClientConnectionError:
            log.error("JWKS fetch failed for %s", self.config.certs_url)
            await _deny(scope, send, 503, "auth temporarily unavailable")
            return
        except PyJWTError as exc:
            log.info("rejected token: %s", type(exc).__name__)
            await _deny(scope, send, 403, "invalid access token")
            return

        scope.setdefault("state", {})["access_claims"] = claims
        await self.app(scope, receive, send)


async def _deny(scope: Scope, send: Send, status: int, detail: str) -> None:
    if scope["type"] == "websocket":
        await send({"type": "websocket.close", "code": 1008})
        return
    body = json.dumps({"detail": detail}).encode()
    await send(
        {
            "type": "http.response.start",
            "status": status,
            "headers": [
                (b"content-type", b"application/json"),
                (b"content-length", str(len(body)).encode()),
                (b"cache-control", b"no-store"),
            ],
        }
    )
    await send({"type": "http.response.body", "body": body})


def install_access_auth(app: Any, **kwargs: Any) -> None:
    """Attach the middleware, or (local only) skip it with ALLOW_INSECURE_DEV=1."""
    if os.environ.get("ALLOW_INSECURE_DEV") == "1":
        if os.environ.get("K_SERVICE"):
            raise RuntimeError("ALLOW_INSECURE_DEV is not permitted on Cloud Run.")
        log.warning("Cloudflare Access auth DISABLED (local dev only).")
        return
    config = AccessConfig.from_env()  # raises before anything is attached
    app.add_middleware(CloudflareAccessMiddleware, config=config, **kwargs)
