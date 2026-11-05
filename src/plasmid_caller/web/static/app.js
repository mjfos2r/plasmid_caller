"use strict";
// Contig IDs and filenames come from uploads: build DOM nodes, never HTML strings.
const view = document.getElementById("view");
let poll = null;

function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (key.startsWith("on")) node.addEventListener(key.slice(2), value);
    else if (value !== null && value !== undefined) node.setAttribute(key, value);
  }
  node.append(...children.filter((child) => child !== null && child !== undefined));
  return node;
}

async function api(path, options) {
  const response = await fetch(path, options);
  const body = await response.json().catch(() => ({}));
  if (!response.ok) throw new Error(body.detail || `Request failed (${response.status})`);
  return body;
}

const number = (value, digits = 0) => value === null || value === undefined ? "" :
  Number(value).toLocaleString(undefined, { minimumFractionDigits: digits, maximumFractionDigits: digits });
const when = (seconds) => new Date(seconds * 1000).toLocaleString();
const statusCell = (job) => el("td", { class: `status-${job.status}` }, job.status);
const showError = (error) => view.replaceChildren(el("p", { class: "error" }, error.message));

function dropZone(onFiles) {
  const input = el("input", { type: "file", multiple: "", accept: ".fa,.fna,.fasta,.gz",
    "aria-label": "Choose assemblies" });
  const zone = el("label", { class: "drop" }, input,
    el("strong", {}, "Drop assemblies here"),
    el("div", { class: "muted" }, "or click to choose .fa, .fna or .fasta files, optionally gzipped"));
  input.addEventListener("change", () => onFiles(input.files));
  zone.addEventListener("dragover", (event) => { event.preventDefault(); zone.classList.add("over"); });
  zone.addEventListener("dragleave", () => zone.classList.remove("over"));
  zone.addEventListener("drop", (event) => {
    event.preventDefault();
    zone.classList.remove("over");
    onFiles(event.dataTransfer.files);
  });
  return zone;
}

async function upload(files, message) {
  if (!files.length) return;
  const form = new FormData();
  for (const file of files) form.append("files", file);
  message.className = "muted";
  message.textContent = `Uploading ${files.length} file(s)…`;
  try {
    const batch = await api("/api/batches", { method: "POST", body: form });
    location.hash = `#/batch/${batch.id}`;
  } catch (error) {
    message.className = "error";
    message.textContent = error.message;
  }
}

async function home() {
  const message = el("p", { "aria-live": "polite" });
  view.replaceChildren(dropZone((files) => upload(files, message)), message);
  const batches = await api("/api/batches");
  if (!batches.length) return;
  const rows = batches.map((batch) => {
    const done = batch.jobs.filter((job) => job.status === "done").length;
    const names = batch.jobs.slice(0, 3).map((job) => job.name).join(", ");
    return el("tr", {},
      el("td", {}, el("a", { href: `#/batch/${batch.id}` }, when(batch.created_at))),
      el("td", {}, names + (batch.jobs.length > 3 ? ` +${batch.jobs.length - 3}` : "")),
      el("td", { class: "num" }, `${done}/${batch.jobs.length}`));
  });
  view.append(el("h2", {}, "Previous batches"), el("div", { class: "scroll" }, el("table", {},
    el("thead", {}, el("tr", {}, el("th", {}, "Submitted"), el("th", {}, "Assemblies"),
      el("th", { class: "num" }, "Done"))),
    el("tbody", {}, ...rows))));
}

async function batchPage(id) {
  const batch = await api(`/api/batches/${id}`);
  const rows = batch.jobs.map((job) => el("tr", {},
    el("td", {}, job.status === "done" ? el("a", { href: `#/job/${job.id}` }, job.name) : job.name),
    statusCell(job),
    el("td", { class: "num" }, number(job.contigs)),
    el("td", { class: "num" }, number(job.bases)),
    el("td", {}, job.error || "")));
  view.replaceChildren(
    el("h1", {}, "Batch"),
    el("p", { class: "muted" }, `Submitted ${when(batch.created_at)} · plasmid_caller ${batch.version}`),
    el("div", { class: "scroll" }, el("table", {},
      el("thead", {}, el("tr", {}, el("th", {}, "Assembly"), el("th", {}, "Status"),
        el("th", { class: "num" }, "Contigs"), el("th", { class: "num" }, "Bases"), el("th", {}, "Note"))),
      el("tbody", {}, ...rows))));
  if (batch.jobs.some((job) => job.status === "queued" || job.status === "running")) {
    poll = setTimeout(() => batchPage(id).catch(showError), 2000);
  }
}

const COLUMNS = [
  { key: "contig_id", label: "Contig" },
  { key: "contig_len", label: "Length", num: true },
  { key: "filtered_call", label: "Call" },
  { key: "genospecies", label: "Genospecies" },
  { key: "plasmid_name_pf32", label: "PF32" },
  { key: "plasmid_name_wp", label: "WP" },
  { key: "strain_wp", label: "WP strain" },
  { key: "overall_percent_identity_wp", label: "WP identity %", num: true, digits: 2 },
  { key: "query_coverage_percent_wp", label: "WP coverage %", num: true, digits: 1 },
];

const isFlagged = (row) => row.genospecies === "ambiguous" || (row.plasmid_name_pf32 &&
  row.plasmid_name_wp && row.plasmid_name_pf32 !== row.plasmid_name_wp);

function contigTable(contigs) {
  let sortKey = "contig_len", ascending = false, query = "";
  const body = el("tbody");
  const render = () => {
    const rows = contigs
      .filter((row) => COLUMNS.some((col) => String(row[col.key] ?? "").toLowerCase().includes(query)))
      .sort((a, b) => {
        const [x, y] = [a[sortKey], b[sortKey]];
        if (x === null) return 1;
        if (y === null) return -1;
        const order = typeof x === "number" ? x - y : String(x).localeCompare(String(y), undefined, { numeric: true });
        return ascending ? order : -order;
      });
    body.replaceChildren(...rows.map((row) => el("tr", { class: isFlagged(row) ? "flagged" : null },
      ...COLUMNS.map((col) => el("td", { class: col.num ? "num" : null },
        col.num ? number(row[col.key], col.digits) : row[col.key] ?? "")))));
  };
  const head = COLUMNS.map((col) => el("th", { class: col.num ? "num" : null }, el("button", {
    onclick: () => {
      ascending = sortKey === col.key ? !ascending : !col.num;
      sortKey = col.key;
      render();
    } }, col.label)));
  const filter = el("input", { class: "filter", type: "search", placeholder: "Filter contigs",
    "aria-label": "Filter contigs", oninput: (event) => { query = event.target.value.toLowerCase(); render(); } });
  render();
  return [filter, el("div", { class: "scroll" }, el("table", {}, el("thead", {}, el("tr", {}, ...head)), body))];
}

async function jobPage(id) {
  const { job, contigs, composition, downloads } = await api(`/api/jobs/${id}/results`);
  const classified = composition.filter((row) => row.genospecies !== "unclassified");
  const summary = classified.map((row) => `${row.genospecies} ${number(row.percent_composition, 1)}%`).join(" · ");
  const flagged = contigs.filter(isFlagged).length;
  view.replaceChildren(
    el("p", {}, el("a", { href: `#/batch/${job.batch_id}` }, "Back to batch")),
    el("h1", {}, job.name),
    el("p", { class: "muted" }, `${number(job.contigs)} contigs · ${number(job.bases)} bp · ${summary || "no classified contigs"}`),
    el("h2", {}, flagged ? `Contigs (${flagged} flagged: PF32/WP disagree or ambiguous genospecies)` : "Contigs"),
    ...contigTable(contigs),
    el("h2", {}, "Downloads"),
    el("div", { class: "links" }, ...downloads.map((name) =>
      el("a", { href: `/api/jobs/${id}/download/${name}` }, name))));
}

function route() {
  clearTimeout(poll);
  const [, page, id] = location.hash.split("/");
  const load = page === "batch" ? batchPage(id) : page === "job" ? jobPage(id) : home();
  load.catch(showError);
}

window.addEventListener("hashchange", route);
route();
