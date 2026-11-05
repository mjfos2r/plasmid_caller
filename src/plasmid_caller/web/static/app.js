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

function compositionBars(jobs) {
  const totals = new Map();
  for (const job of jobs) for (const row of job.composition) {
    if (row.percent_composition > 0) totals.set(row.genospecies, (totals.get(row.genospecies) ?? 0) + row.percent_composition);
  }
  if (!totals.size) return [];
  // Colour follows the species for the whole batch; anything past the seventh folds into Other.
  const named = [...totals.keys()].sort((a, b) => totals.get(b) - totals.get(a) || a.localeCompare(b)).slice(0, 7);
  const slot = (name) => named.includes(name) ? `var(--c${named.indexOf(name)})` : "var(--muted)";
  const swatch = (label, colour) => {
    const mark = el("span", { class: "swatch" });
    mark.style.background = colour;
    return el("span", {}, mark, label);
  };
  const legend = el("div", { class: "legend" }, ...named.map((name) => swatch(name, slot(name))),
    totals.size > named.length ? swatch("Other", "var(--muted)") : null);
  const rows = jobs.map((job) => {
    const parts = job.composition.filter((row) => row.percent_composition > 0)
      .sort((a, b) => b.percent_composition - a.percent_composition);
    const text = parts.map((row) => `${row.genospecies} ${number(row.percent_composition, 1)}%`).join(", ");
    const stack = el("div", { class: "stack", role: "img", "aria-label": text || "No classified contigs" },
      ...parts.map((row) => {
        const part = el("span", { title: `${row.genospecies}: ${number(row.percent_composition, 1)}%` });
        part.style.flexGrow = row.percent_composition;
        part.style.background = slot(row.genospecies);
        return part;
      }));
    return el("div", { class: "track composition" }, el("div", { class: "track-label", title: job.name },
      el("a", { href: `#/job/${job.id}` }, job.name)),
      el("div", {}, stack, el("div", { class: "muted" }, text || "No classified contigs")));
  });
  return [el("h2", {}, "Genospecies composition"), legend, ...rows];
}

function repliconTable(jobs) {
  const names = [...new Set(jobs.flatMap((job) => Object.keys(job.replicons)))].sort(chromosomeFirst);
  if (!names.length) return [];
  const peak = Object.fromEntries(names.map((name) => [name, Math.max(...jobs.map((job) => job.replicons[name]?.bp ?? 0))]));
  const cell = (job, name) => {
    const entry = job.replicons[name];
    if (!entry) return el("td");
    const node = el("td", { class: "num", title: `${name}: ${number(entry.bp)} bp in ${number(entry.contigs)} contig(s)` },
      number(entry.bp / 1000, 1));
    node.style.setProperty("--shade", `${10 + entry.bp / peak[name] * 35}%`);
    return node;
  };
  return [el("h2", {}, "Called replicons (kb)"),
    el("p", { class: "muted" }, "Bases called to each replicon, shaded within each column. A multi-locus contig counts towards each of its names."),
    el("div", { class: "scroll matrix" }, el("table", {},
      el("thead", {}, el("tr", {}, el("th", { scope: "col" }, "Assembly"),
        ...names.map((name) => el("th", { class: "num" }, name)))),
      el("tbody", {}, ...jobs.map((job) => el("tr", {},
        el("th", { scope: "row", title: job.name }, el("a", { href: `#/job/${job.id}` }, job.name)),
        ...names.map((name) => cell(job, name)))))))];
}

async function batchPage(id) {
  const [batch, overview] = await Promise.all([api(`/api/batches/${id}`), api(`/api/batches/${id}/overview`)]);
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
      el("tbody", {}, ...rows))),
    ...compositionBars(overview), ...repliconTable(overview));
  if (batch.jobs.some((job) => job.status === "queued" || job.status === "running")) {
    poll = setTimeout(() => batchPage(id).catch(showError), 2000);
  }
}

const contigHref = (jobId, contigId) => `#/job/${jobId}/contig/${encodeURIComponent(contigId)}`;

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
  { key: "fusion_signal", label: "Fusion candidate" },
];

const isFlagged = (row) => row.genospecies === "ambiguous" || (row.plasmid_name_pf32 &&
  row.plasmid_name_wp && row.plasmid_name_pf32 !== row.plasmid_name_wp);

function contigTable(contigs, jobId) {
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
        col.key === "contig_id" ? el("a", { href: contigHref(jobId, row.contig_id) }, row.contig_id) :
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

const chromosomeFirst = (a, b) => (b === "chromosome") - (a === "chromosome") ||
  a.localeCompare(b, undefined, { numeric: true });
// Callable hits rank above grey ones, so a short uncallable 100% can't lead the table.
const columnRank = (cell) => !cell ? -1 : cell.overall_percent_identity + (cell.passes ? 1000 : 0);
const identityShade = (identity) => `${Math.max(0, Math.min(1, (identity - 80) / 20)) * 45}%`;

function identityCell(label, row, otherLabel) {
  const cell = row.get(label), fused = row.get("fusion");
  const isPart = Boolean(fused?.passes && fused.parts.includes(label));
  if (!cell) return isPart ? el("td", { class: "num part", title: `Part of ${fused.fusion}, with no unfused hit` }, "–") : el("td");
  const lines = [
    `${cell.label}: ${number(cell.overall_percent_identity, 2)}% identity`,
    `Covers ${number(cell.query_coverage_percent, 1)}% of the contig (${number(cell.query_covered_length)} bp)`,
    `${otherLabel}: ${cell.other ?? "unknown"}`,
    `Reference: ${cell.strain ?? "unknown strain"} ${cell.plasmid_id ?? ""}`,
  ];
  if (cell.fusion) {
    lines[0] = `${cell.fusion}: ${number(cell.overall_percent_identity, 2)}% identity`;
    lines.push(`Parts: ${cell.parts.map((part) => {
      const hit = row.get(part);
      return `${part} ${hit ? `${number(hit.overall_percent_identity, 2)}%` : "no unfused hit"}`;
    }).join(", ")}`);
  }
  if (isPart) lines.push(`Part of ${fused.fusion}`);
  if (cell.top) lines.push("Highest-scoring WP hit on this contig");
  if (!cell.passes) lines.push("Below the coverage gates: cannot be called");
  const node = el("td", { class: `num${cell.top ? " top" : ""}${cell.passes ? "" : " weak"}${isPart ? " part" : ""}`,
    title: lines.join("\n") }, number(cell.overall_percent_identity, 2),
    cell.fusion ? el("div", { class: "muted" }, cell.fusion) : null);
  if (cell.passes) node.style.setProperty("--shade", identityShade(cell.overall_percent_identity));
  return node;
}

function identityTable(contigs, cells, otherLabel, jobId) {
  if (!cells.length) return [el("p", { class: "muted" }, "No WP hits for this assembly.")];
  const labels = [...new Set(cells.flatMap((cell) => [cell.label, ...cell.parts]))].sort(chromosomeFirst);
  const hits = new Map(contigs.map((row) => [row.contig_id, new Map()]));
  for (const cell of cells) hits.get(cell.contig_id)?.set(cell.label, cell);
  const identity = (row, label) => hits.get(row.contig_id).get(label)?.overall_percent_identity ?? -1;
  const peak = new Map(labels.map((label) => [label, -1]));
  for (const cell of cells) peak.set(cell.label, Math.max(peak.get(cell.label), columnRank(cell)));
  let sortLabel = null, query = "", anchor = null;
  const body = el("tbody"), headRow = el("tr");
  const anchorOption = el("option", { value: "contig" });
  const order = el("select", { "aria-label": "Column order", onchange: () => {
    if (order.value !== "contig") setAnchor(null);
    render();
  } }, el("option", { value: "identity" }, "Highest identity"), el("option", { value: "name" }, "Alphabetical"));
  const setAnchor = (contigId) => {
    anchor = contigId;
    if (anchor === null) {
      if (order.value === "contig") order.value = "identity";
      anchorOption.remove();
    } else {
      anchorOption.textContent = `Contig ${anchor}`;
      order.append(anchorOption);
      order.value = "contig";
    }
  };
  // The fusion column stays first in every order, beside the contig names.
  const pinned = labels.filter((label) => label === "fusion");
  const others = labels.filter((label) => label !== "fusion");
  const columns = () => {
    if (order.value === "name") return [...pinned, ...others];
    const rank = (label) => anchor === null ? peak.get(label) : columnRank(hits.get(anchor).get(label));
    return [...pinned, ...[...others].sort((a, b) => rank(b) - rank(a) || chromosomeFirst(a, b))];
  };
  const pin = (label, node) => {
    if (label === "fusion") node.classList.add("pin");
    return node;
  };
  const render = () => {
    const shown = columns();
    headRow.replaceChildren(el("th", { scope: "col" }, "Contig"), ...shown.map((label) =>
      pin(label, el("th", { class: "num" }, el("button", {
        onclick: () => { sortLabel = sortLabel === label ? null : label; render(); } }, label)))));
    const rows = contigs
      .filter((row) => `${row.contig_id} ${row.filtered_call ?? ""}`.toLowerCase().includes(query))
      .sort((a, b) => sortLabel === null ? b.contig_len - a.contig_len : identity(b, sortLabel) - identity(a, sortLabel));
    body.replaceChildren(...rows.map((row) => el("tr", { class: row.contig_id === anchor ? "anchor" : null },
      el("th", { scope: "row", title: row.contig_id },
        el("button", { "aria-pressed": row.contig_id === anchor, onclick: () => {
          setAnchor(row.contig_id === anchor ? null : row.contig_id);
          render();
        } }, row.contig_id),
        el("div", { class: "muted" }, `${row.filtered_call ?? "no call"} · ${number(row.contig_len)} bp · `,
          el("a", { href: contigHref(jobId, row.contig_id) }, "open"))),
      ...shown.map((label) => pin(label, identityCell(label, hits.get(row.contig_id), otherLabel))))));
  };
  const filter = el("input", { class: "filter", type: "search", placeholder: "Filter contigs",
    "aria-label": "Filter contigs", oninput: (event) => { query = event.target.value.toLowerCase(); render(); } });
  render();
  return [el("div", { class: "controls" }, filter, el("label", {}, "Columns ", order)),
    el("div", { class: "scroll matrix" }, el("table", {}, el("thead", {}, headRow), body))];
}

function tabs(panels) {
  const bodies = panels.map(([, nodes]) => el("div", {}, ...nodes));
  const buttons = panels.map(([label], index) => el("button", { class: "tab", onclick: () => show(index) }, label));
  const show = (index) => {
    bodies.forEach((body, i) => { body.hidden = i !== index; });
    buttons.forEach((button, i) => button.setAttribute("aria-pressed", i === index));
  };
  show(0);
  return [el("div", { class: "tabs" }, ...buttons), ...bodies];
}

async function jobPage(id) {
  const [{ job, contigs, composition, downloads }, identity] = await Promise.all([
    api(`/api/jobs/${id}/results`), api(`/api/jobs/${id}/identity`)]);
  for (const row of contigs) {
    const signals = [row.multiple_loci_pf32 && "PF32", identity.candidates[row.contig_id] && "WP"].filter(Boolean);
    row.fusion_signal = signals.length ? signals.join(" + ") : null;
  }
  const classified = composition.filter((row) => row.genospecies !== "unclassified");
  const summary = classified.map((row) => `${row.genospecies} ${number(row.percent_composition, 1)}%`).join(" · ");
  const flagged = contigs.filter(isFlagged).length;
  view.replaceChildren(
    el("p", {}, el("a", { href: `#/batch/${job.batch_id}` }, "Back to batch")),
    el("h1", {}, job.name),
    el("p", { class: "muted" }, `${number(job.contigs)} contigs · ${number(job.bases)} bp · ${summary || "no classified contigs"}`),
    ...tabs([
      ["Contigs", [
        el("h2", {}, flagged ? `Contigs (${flagged} flagged: PF32/WP disagree or ambiguous genospecies)` : "Contigs"),
        ...contigTable(contigs, id)]],
      ["Plasmid identity", identityTable(contigs, identity.plasmid, "Genospecies", id)],
      ["Genospecies identity", identityTable(contigs, identity.genospecies, "Plasmid", id)],
    ]),
    ...(job.provenance ? [el("p", { class: "muted" }, `plasmid_caller ${job.provenance.caller_version} · BLAST ${job.provenance.blast_version} · ` +
      `database ${job.provenance.database_sha256.slice(0, 12)} · run ${when(job.finished_at)}`)] : []),
    el("h2", {}, "Downloads"),
    el("div", { class: "links" }, ...downloads.map((name) =>
      el("a", { href: `/api/jobs/${id}/download/${name}` }, name))));
}

function hitTrack(label, hit, length, { detail, weak = false }) {
  const bar = el("div", { class: "bar" });
  for (const [start, end] of hit.query_intervals) {
    const span = el("span", { class: weak ? "weak" : null });
    span.style.left = `${(start - 1) / length * 100}%`;
    span.style.width = `${Math.max((end - start + 1) / length * 100, 0.2)}%`;
    if (!weak) span.style.setProperty("--shade", `${30 + parseFloat(identityShade(hit.overall_percent_identity)) * 1.5}%`);
    bar.append(span);
  }
  return el("div", { class: "track", title: detail.join("\n") },
    el("div", { class: "track-label" }, label), bar,
    el("div", { class: "num muted" }, `${number(hit.overall_percent_identity, 2)}%`));
}

function wpTracks(cells, length, name) {
  if (!cells.length) return [el("p", { class: "muted" }, "No WP hits on this contig.")];
  const ranked = [...cells].sort((a, b) => b.passes - a.passes || b.score - a.score);
  const list = el("div");
  const render = (showGrey) => list.replaceChildren(...ranked.filter((cell) => showGrey || cell.passes).map((cell) =>
    hitTrack(name(cell), cell, length, { weak: !cell.passes, detail: [
      `${number(cell.overall_percent_identity, 2)}% identity over ${number(cell.query_coverage_percent, 1)}% of the contig (${number(cell.query_covered_length)} bp)`,
      `Reference: ${cell.strain ?? "unknown strain"} ${cell.plasmid_id ?? ""}, ${number(cell.ref_covered_length)} of ${number(cell.ref_length)} bp covered`,
      ...(cell.fusion ? [`Fused reference: ${cell.parts.join(" + ")}`] : []),
      ...(cell.passes ? [] : ["Below the coverage gates: cannot be called"]),
    ] })));
  render(false);
  const grey = cells.length - cells.filter((cell) => cell.passes).length;
  const toggle = grey ? el("label", { class: "muted" }, el("input", { type: "checkbox",
    onchange: (event) => render(event.target.checked) }), ` Show ${grey} hits below the coverage gates`) : null;
  return [list, toggle];
}

function fusionPanel(contig, loci, tiling) {
  const lines = [];
  if (contig.multiple_loci_pf32) {
    lines.push(`PF32 loci: ${loci.map((locus) =>
      `${locus.plasmid_name} at ${number(locus.locus_start)}–${number(locus.locus_end)}`).join(", ")}`);
  }
  if (tiling.length >= 2) {
    lines.push(`WP tiling: ${tiling.map(([name, bp], index) =>
      `${name} (${index ? "+" : ""}${number(bp)} bp)`).join(", ")}`);
  }
  if (!lines.length) return el("p", { class: "muted" }, "No fusion signals on this contig.");
  return el("div", { class: "notice" }, el("strong", {}, "Fusion candidate: needs validation"),
    ...lines.map((line) => el("div", {}, line)));
}

async function contigPage(jobId, contigId) {
  const { job, contig, rule, plasmids, genospecies, pf32, loci, tiling } =
    await api(`/api/jobs/${jobId}/contig?contig_id=${encodeURIComponent(contigId)}`);
  const length = contig.contig_len;
  const species = genospecies.filter((cell) => cell.passes).sort((a, b) => b.score - a.score);
  const margin = species.length > 1 ? `${species[0].label} leads ${species[1].label} on coverage × identity: ` +
    `${number(species[0].score / 100, 1)} against ${number(species[1].score / 100, 1)}.` : null;
  const pf32Tracks = pf32.map((hit) => hitTrack(`${hit.plasmid_name} · ${hit.strain ?? "?"}`, hit, length, { detail: [
    `${number(hit.overall_percent_identity, 2)}% identity over ${number(hit.query_covered_length)} bp`,
    `Reference: ${hit.strain ?? "unknown strain"} ${hit.plasmid_id ?? ""}`] }));
  view.replaceChildren(
    el("p", {}, el("a", { href: `#/job/${jobId}` }, `Back to ${job.name}`)),
    el("h1", {}, contig.contig_id),
    el("p", { class: "muted" }, `${number(length)} bp · ${contig.genospecies ?? "no genospecies"}`),
    el("p", {}, "Call: ", el("strong", {}, contig.filtered_call ?? "none"), rule ? ` · ${rule}` : ""),
    fusionPanel(contig, loci, tiling),
    el("div", { class: "axis muted" }, el("span", {}, "1"), el("span", {}, `${number(length)} bp`)),
    el("h2", {}, "PF32 hits"),
    ...(pf32Tracks.length ? pf32Tracks : [el("p", { class: "muted" }, "No PF32 hits on this contig.")]),
    ...tabs([
      ["WP by plasmid", wpTracks(plasmids, length, (cell) => cell.label)],
      ["WP by genospecies", [margin ? el("p", { class: "muted" }, margin) : null,
        ...wpTracks(genospecies, length, (cell) => cell.label)].filter(Boolean)],
    ]));
}

function route() {
  clearTimeout(poll);
  const [, page, id, sub, contig] = location.hash.split("/");
  const load = page === "batch" ? batchPage(id) :
    page === "job" && sub === "contig" ? contigPage(id, decodeURIComponent(contig)) :
    page === "job" ? jobPage(id) : home();
  load.catch(showError);
}

window.addEventListener("hashchange", route);
route();
