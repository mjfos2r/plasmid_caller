#!/usr/bin/env bash
# Build the web portal via Cloud Build (no local Docker needed) and deploy it to Cloud Run.
#
# Usage:
#   PROJECT_ID=... BUCKET_NAME=... CF_TEAM_DOMAIN=... CF_ACCESS_AUD=... PRIVATE_HOSTNAME=... \
#     deploy/build_and_deploy.sh <tag>      e.g. v0.4
#
# Nothing project-specific is hardcoded here on purpose -- CF_ACCESS_AUD especially is sensitive
# and must never be committed. Export these in your own shell (or source a local, gitignored
# env file) before running. See deploy/README.md for what each one means.
set -Eeuo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

: "${PROJECT_ID:?Set PROJECT_ID (e.g. mjf-personal)}"
: "${BUCKET_NAME:?Set BUCKET_NAME, the GCS bucket for called genomes}"
: "${CF_TEAM_DOMAIN:?Set CF_TEAM_DOMAIN (e.g. yourteam.cloudflareaccess.com)}"
: "${CF_ACCESS_AUD:?Set CF_ACCESS_AUD, the Access application's Audience tag}"
: "${PRIVATE_HOSTNAME:?Set PRIVATE_HOSTNAME, the SSO-gated hostname (Cloud Run's own PUBLIC_HOSTNAME env var)}"

REGION=us-east1
SERVICE=plasmid-portal
IMAGE_REPO=portal/web
RUNTIME_SA="portal-runtime@${PROJECT_ID}.iam.gserviceaccount.com"

tag="${1:?usage: $0 <tag>, e.g. v0.4}"
img="${REGION}-docker.pkg.dev/${PROJECT_ID}/${IMAGE_REPO}:${tag}"
revision="$(git rev-parse HEAD)"

echo "==> Building ${img} from ${revision}"
gcloud builds submit --config=deploy/cloudbuild.web.yaml \
  --substitutions=_IMG="${img}",_REVISION="${revision}",_WEB_VERSION="${tag}" .

echo "==> Deploying ${img} to ${SERVICE}"
gcloud run deploy "${SERVICE}" \
  --image="${img}" --region="${REGION}" \
  --service-account="${RUNTIME_SA}" \
  --allow-unauthenticated \
  --max-instances=1 --concurrency=10 --no-cpu-throttling \
  --cpu=2 --memory=8Gi --timeout=300 \
  --set-env-vars="CF_TEAM_DOMAIN=${CF_TEAM_DOMAIN},CF_ACCESS_AUD=${CF_ACCESS_AUD},PUBLIC_HOSTNAME=${PRIVATE_HOSTNAME},BUCKET_NAME=${BUCKET_NAME}"

echo "==> Done: ${img}"
