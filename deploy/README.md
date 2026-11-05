# Deploying the web portal

Browser → Cloudflare Access → Worker → Cloud Run. The `run.app` URL is public, so the
JWT check in the app (`web/auth.py`) is the real lock. Do not share the link until every
check below passes.

Replace `<PROJECT_ID>`, `<TEAM>`, `<HOSTNAME>`, `<AUD>` and `<TAG>`.

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
  --build-arg SOURCE_REVISION=$(git rev-parse HEAD) -t $IMG .
docker push $IMG
```

`--platform linux/amd64` is required on Apple Silicon, because Cloud Run needs an amd64 image.
Emulated builds are slower.

## 3. Deploy

```bash
gcloud run deploy plasmid-portal \
  --image=$IMG --region=us-east1 \
  --service-account=portal-runtime@<PROJECT_ID>.iam.gserviceaccount.com \
  --allow-unauthenticated \
  --max-instances=1 --concurrency=10 --no-cpu-throttling \
  --cpu=2 --memory=4Gi --timeout=300 \
  --set-env-vars=CF_TEAM_DOMAIN=<TEAM>.cloudflareaccess.com,CF_ACCESS_AUD=<AUD>,PUBLIC_HOSTNAME=<HOSTNAME>
```

These differ from a typical stateless service on purpose:

- `--max-instances=1`: jobs live in one process, so a second instance would answer polls with 404.
- `--concurrency=10`: pages keep polling while a job runs.
- `--no-cpu-throttling`: jobs run after the upload response, when Cloud Run would otherwise stop the CPU.

Results are held in the instance's memory and are lost when it stops or on any deploy.
Tell users to download what they need.

## 4. Worker

Create a Worker from `worker.js`, set the variable `ORIGIN` to the `run.app` URL, and add
`<HOSTNAME>` as its custom domain. It must match the Access application's hostname.

## 5. Check before sharing

Lock:

```bash
curl -i https://<service>.run.app/                                                    # 403
curl -i -H 'Cf-Access-Jwt-Assertion: a.b.c' https://<service>.run.app/                # 403
curl -i https://<service>.run.app/api/version                                         # 403
```

- [ ] Private window to `https://<HOSTNAME>` shows the Cloudflare login
- [ ] An email outside the allowlist is denied
- [ ] An allowlisted email loads the app, and `/api/version` shows the expected stamps
- [ ] An assembly uploads through `<HOSTNAME>` and finishes

Whether Access passes `Cf-Access-Jwt-Assertion` to a Worker custom domain is unconfirmed.
If the allowlisted login gets a 403 from the app, start there.

Parity, against the deployed portal. Run this on a Linux machine with the CLI installed, since it
compares with a local run, and use a real token from the browser's request headers:

```bash
PORTAL_URL=https://<HOSTNAME> PORTAL_HEADERS='{"Cf-Access-Jwt-Assertion": "<JWT>"}' \
PARITY_FASTAS="a.fna b.fna" uv run --all-extras pytest tests/test_parity.py -s
```

## Update or roll back

Push a new tag and rerun step 3 with it, or redeploy the previous tag. Either one clears all results.
