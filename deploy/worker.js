// Thin proxy, shared by both the private and public custom domains. Set ORIGIN to the
// https://...run.app URL, and PUBLIC_SITE_HOSTNAME to whichever domain should serve the public
// SPA at "/" (e.g. plasmids.syntenous.cc) -- distinct from Cloud Run's own PUBLIC_HOSTNAME env
// var, which is the *private* hostname, used there for a different purpose (Origin validation).
//
// The app's /public route exists either way; this just means visitors never see "/public" in
// the address bar on that domain. Nothing else needs rewriting -- the public SPA routes
// entirely via the URL fragment (#/genome/...), which browsers never send to the server, so
// every reload still hits this same "/" rewrite regardless of which in-app page was open.
export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    if (env.PUBLIC_SITE_HOSTNAME && url.hostname === env.PUBLIC_SITE_HOSTNAME && url.pathname === "/") {
      url.pathname = "/public";
    }
    const target = new URL(url.pathname + url.search, env.ORIGIN);
    return fetch(new Request(target, request));
  },
};
