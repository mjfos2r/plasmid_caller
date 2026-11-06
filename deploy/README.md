# Deploying the web portal

Two hostnames, one Worker, one Cloud Run service:

- `<PRIVATE_HOSTNAME>` -- the authenticated app (upload, manage visibility, full results).
  Browser → Cloudflare Access → Worker → Cloud Run.
- `<PUBLIC_HOSTNAME>` -- the read-only `/public` surface, with no Cloudflare Access in front of it
  at all. Browser → Worker → Cloud Run.

The `run.app` URL itself is public either way, so the JWT check in the app (`web/auth.py`) is the
real lock, not Cloudflare Access. Do not share either link until every check in step 5 passes.

Replace `<PROJECT_ID>`, `<TEAM>`, `<PRIVATE_HOSTNAME>`, `<PUBLIC_HOSTNAME>`, `<AUD>`, `<BUCKET>`
and `<TAG>`.

## 1. Get the code and use the right account

Needs Docker with BuildKit (Docker Desktop has it) and the gcloud CLI. The BLAST binaries and
seqkit are fetched during the build, and the reference database is tracked in git.

```bash
git clone https://github.com/mjfos2r/plasmid_caller.git
cd plasmid_caller
git switch web-viewer
```

```bash
gcloud auth login mikejoshfoster@gmail.com
gcloud config set project <PROJECT_ID>
gcloud config list            # confirm before continuing
gcloud auth configure-docker us-east1-docker.pkg.dev
```

## 2. Build and push

```bash
IMG=us-east1-docker.pkg.dev/<PROJECT_ID>/portal/web:<TAG>
docker build --platform linux/amd64 -f Dockerfile.web \
  --build-arg SOURCE_REVISION=$(git rev-parse HEAD) --build-arg WEB_VERSION=<TAG> -t $IMG .
docker push $IMG
```

`--platform linux/amd64` is required on Apple Silicon, because Cloud Run needs an amd64 image.
Emulated builds are slower. Building via Cloud Build instead of local Docker sidesteps this --
submit with a config that builds `-f Dockerfile.web` explicitly, since the repo's default
`Dockerfile` is for the CLI image, not the portal, and neither `gcloud run deploy --source` nor
the `gcloud builds submit --tag` shorthand support a custom Dockerfile path.

## 3. Deploy

```bash
gcloud run deploy plasmid-portal \
  --image=$IMG --region=us-east1 \
  --service-account=portal-runtime@<PROJECT_ID>.iam.gserviceaccount.com \
  --allow-unauthenticated \
  --max-instances=1 --concurrency=10 --no-cpu-throttling \
  --cpu=2 --memory=8Gi --timeout=300 \
  --set-env-vars=CF_TEAM_DOMAIN=<TEAM>.cloudflareaccess.com,CF_ACCESS_AUD=<AUD>,PUBLIC_HOSTNAME=<PRIVATE_HOSTNAME>,BUCKET_NAME=<BUCKET>
```

`PUBLIC_HOSTNAME` must be the *private* hostname, not the public one: behind the Worker, `Host`
always shows up as the `run.app` address, so the app validates uploads' `Origin` header against
`PUBLIC_HOSTNAME` instead -- and uploads only ever happen through the authenticated side.

These differ from a typical stateless service on purpose:

- `--max-instances=1`: jobs live in one process, so a second instance would answer polls with 404.
- `--concurrency=10`: pages keep polling while a job runs.
- `--no-cpu-throttling`: jobs run after the upload response, when Cloud Run would otherwise stop the CPU.

Batch/job history is still held in the instance's memory and is lost when it stops or on any
deploy -- tell users to download what they need. Called genomes and the dedup index are not:
they're written through to `<BUCKET>` (give `portal-runtime@<PROJECT_ID>.iam.gserviceaccount.com`
`roles/storage.objectAdmin` on it), so a redeploy no longer loses them, and an already-called
genome short-circuits to a `found` result instead of rerunning BLAST.

## 4. Cloudflare: two hostnames, one Worker

Create a Worker from `worker.js` and set the variable `ORIGIN` to the `run.app` URL. It's almost a
plain passthrough -- it doesn't care which hostname a request came in on for auth purposes, so one
Worker serves both:

- Add `<PRIVATE_HOSTNAME>` as a custom domain on the Worker, and put a Cloudflare Access
  application in front of it with your login policy. This must match the Access application's
  hostname.
- Add `<PUBLIC_HOSTNAME>` as a second custom domain on the *same* Worker, with **no** Access
  application in front of it. Also set the Worker variable `PUBLIC_SITE_HOSTNAME` to this same
  hostname -- **not** the same thing as Cloud Run's `PUBLIC_HOSTNAME` env var from step 3, which
  is the *private* hostname. The Worker uses it to rewrite `/` to `/public` for this domain only,
  so visitors never see `/public` in the address bar; everything past that is client-side hash
  routing (`#/genome/...`), which the Worker never sees, so nothing else needs rewriting.

The split needs no app code to enforce it: Cloudflare Access only attaches a
`Cf-Access-Jwt-Assertion` header to requests on a hostname that has an Access application, and the
app's own middleware already only treats `/public` and `/api/public/*` as not requiring one. So
`<PUBLIC_HOSTNAME>` can only ever reach those paths -- everything else 403s with no JWT present,
the same way a request straight to the bare `run.app` URL would.

## 5. Check before sharing

Lock, on the private hostname:

```bash
curl -i https://<service>.run.app/                                                    # 403
curl -i -H 'Cf-Access-Jwt-Assertion: a.b.c' https://<service>.run.app/                # 403
curl -i https://<service>.run.app/api/version                                         # 403
```

- [ ] Private window to `https://<PRIVATE_HOSTNAME>` shows the Cloudflare login
- [ ] An email outside the allowlist is denied
- [ ] An allowlisted email loads the app, and `/api/version` shows the expected stamps
- [ ] An assembly uploads through `<PRIVATE_HOSTNAME>` and finishes

Open, on the public hostname:

- [ ] `https://<PUBLIC_HOSTNAME>/public` loads with **no** login prompt
- [ ] `https://<PUBLIC_HOSTNAME>/` and `https://<PUBLIC_HOSTNAME>/api/batches` are blocked (403),
      also with no login prompt -- confirms it genuinely can't reach anything private, not just
      that the UI doesn't link there
- [ ] A genome/batch toggled visible from `<PRIVATE_HOSTNAME>` appears at
      `https://<PUBLIC_HOSTNAME>/public`

Whether Access passes `Cf-Access-Jwt-Assertion` to a Worker custom domain is unconfirmed.
If the allowlisted login gets a 403 from the app, start there.

Parity, against the deployed portal. Run this on a Linux machine with the CLI installed, since it
compares with a local run, and use a real token from the browser's request headers:

```bash
PORTAL_URL=https://<PRIVATE_HOSTNAME> PORTAL_HEADERS='{"Cf-Access-Jwt-Assertion": "<JWT>"}' \
PARITY_FASTAS="a.fna b.fna" uv run --all-extras pytest tests/test_parity.py -s
```

## Update or roll back

Push a new tag and rerun step 3 with it, or redeploy the previous tag. Either one clears batch/job
history, but not called genomes -- those live in `<BUCKET>`, independent of the running revision.
