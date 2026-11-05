// Thin proxy from the public hostname to Cloud Run. Set ORIGIN to the https://...run.app URL.
export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    const target = new URL(url.pathname + url.search, env.ORIGIN);
    return fetch(new Request(target, request));
  },
};
