"use strict";
// Shared between the authenticated SPA (app.js) and the public read-only one (public.js).
// Contig IDs and filenames come from uploads: build DOM nodes, never HTML strings.

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
const showError = (view, error) => view.replaceChildren(el("p", { class: "error" }, error.message));

// Real browser history, not a hardcoded destination -- so it's correct from wherever the user
// actually came from (a deep link, a different tab, hash-routing within the SPA), and every
// page gets one for free instead of each page hand-computing its own "back to X" target.
function backButton() {
  return el("p", {}, el("a", { href: "#", onclick: (event) => { event.preventDefault(); history.back(); } }, "← Back"));
}

function sortableTable(rows, columns, initialSort) {
  let sortKey = initialSort, ascending = true, query = "";
  const body = el("tbody");
  const render = () => {
    const shown = rows
      .filter((row) => columns.some((col) => String(row[col.key] ?? "").toLowerCase().includes(query)))
      .sort((a, b) => {
        const [x, y] = [a[sortKey], b[sortKey]];
        const order = typeof x === "number" ? x - y :
          String(x ?? "").localeCompare(String(y ?? ""), undefined, { numeric: true });
        return ascending ? order : -order;
      });
    body.replaceChildren(...shown.map((row) => el("tr", {},
      ...columns.map((col) => el("td", { class: col.num ? "num" : null },
        col.render ? col.render(row) : col.num ? number(row[col.key]) : row[col.key] ?? "")))));
  };
  const head = columns.map((col) => el("th", { class: col.num ? "num" : null }, el("button", {
    onclick: () => {
      ascending = sortKey === col.key ? !ascending : true;
      sortKey = col.key;
      render();
    } }, col.label)));
  const filter = el("input", { class: "filter", type: "search", placeholder: "Filter",
    "aria-label": "Filter", oninput: (event) => { query = event.target.value.toLowerCase(); render(); } });
  render();
  return [filter, el("div", { class: "scroll" }, el("table", {}, el("thead", {}, el("tr", {}, ...head)), body))];
}

// `focus` ("plasmid" or "genospecies") tells the contig page which WP tab to open on, so a link
// from the genospecies identity table doesn't land on the plasmid tracks by default.
const contigHref = (base, jobId, contigId, focus) =>
  `${base}/${jobId}/contig/${encodeURIComponent(contigId)}${focus ? `/${focus}` : ""}`;
const labelHref = (base, jobId, by, label) => `${base}/${jobId}/${by}/${encodeURIComponent(label)}`;

const ncbiNuccoreUrl = (accession) => `https://www.ncbi.nlm.nih.gov/nuccore/${encodeURIComponent(accession)}`;
const ncbiAssemblyUrl = (accession) => `https://www.ncbi.nlm.nih.gov/datasets/genome/${encodeURIComponent(accession)}/`;

function copyFastaButton(url) {
  const button = el("button", { class: "upload-toggle", onclick: async () => {
    try {
      const response = await fetch(url);
      if (!response.ok) throw new Error(`Request failed (${response.status})`);
      await navigator.clipboard.writeText(await response.text());
      button.textContent = "Copied!";
      setTimeout(() => { button.textContent = "Copy FASTA"; }, 1500);
    } catch (error) {
      alert(error.message);
    }
  } }, "Copy FASTA");
  return button;
}

// All contigs that hit one plasmid, or were called to one genospecies -- reached by clicking
// that name on a single contig's track view. Reuses the job/genome's own identity cells, just
// filtered to one label; no new API shape.
function labelContigsTable(cells, by, contigLengths, contigHrefFor) {
  // Same track visualization as a single contig's own hit list (hitTrack/wpTracks), just one
  // row per contig instead of one row per hit -- each still scaled to that contig's own length.
  const otherLabel = by === "plasmid" ? "Genospecies" : "Plasmid";
  const ranked = [...cells].sort((a, b) => b.score - a.score);
  const tracks = ranked.map((cell) => hitTrack(
    el("a", { href: contigHrefFor(cell.contig_id) }, cell.contig_id),
    cell, contigLengths.get(cell.contig_id) ?? cell.query_covered_length,
    { weak: !cell.passes, detail: [
      `${number(cell.overall_percent_identity, 2)}% identity over ${number(cell.query_coverage_percent, 1)}% of the contig (${number(cell.query_covered_length)} bp)`,
      `${otherLabel}: ${cell.other ?? "unknown"}`,
      `Reference: ${cell.strain ?? "unknown strain"} ${cell.plasmid_id ?? ""}`,
    ] }));
  return [el("p", { class: "muted" }, `${number(cells.length)} contig(s).`), ...tracks];
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
  { key: "fusion_signal", label: "Fusion candidate" },
];

const isFlagged = (row) => row.genospecies === "ambiguous" || (row.plasmid_name_pf32 &&
  row.plasmid_name_wp && row.plasmid_name_pf32 !== row.plasmid_name_wp);

function contigTable(contigs, contigHrefFor) {
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
        col.key === "contig_id" ? el("a", { href: contigHrefFor(row.contig_id) }, row.contig_id) :
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
const coverageRank = (cell) => !cell ? -1 : cell.query_coverage_percent + (cell.passes ? 1000 : 0);
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
    el("div", { class: "muted" }, `${number(cell.query_coverage_percent, 1)}% · ${number(cell.query_covered_length)}bp`),
    cell.fusion ? el("div", { class: "muted" }, cell.fusion) : null);
  if (cell.passes) node.style.setProperty("--shade", identityShade(cell.overall_percent_identity));
  return node;
}

function identityTable(contigs, cells, otherLabel, by, contigHrefFor) {
  if (!cells.length) return [el("p", { class: "muted" }, "No WP hits for this assembly.")];
  const labels = [...new Set(cells.flatMap((cell) => [cell.label, ...cell.parts]))].sort(chromosomeFirst);
  const hits = new Map(contigs.map((row) => [row.contig_id, new Map()]));
  for (const cell of cells) hits.get(cell.contig_id)?.set(cell.label, cell);
  const identity = (row, label) => hits.get(row.contig_id).get(label)?.overall_percent_identity ?? -1;
  const peak = new Map(labels.map((label) => [label, -1]));
  const peakCoverage = new Map(labels.map((label) => [label, -1]));
  for (const cell of cells) {
    peak.set(cell.label, Math.max(peak.get(cell.label), columnRank(cell)));
    peakCoverage.set(cell.label, Math.max(peakCoverage.get(cell.label), coverageRank(cell)));
  }
  let sortLabel = null, query = "", anchor = null;
  const body = el("tbody"), headRow = el("tr");
  const anchorOption = el("option", { value: "contig" });
  const order = el("select", { "aria-label": "Column order", onchange: () => {
    if (order.value !== "contig") setAnchor(null);
    render();
  } }, el("option", { value: "identity" }, "Highest identity"), el("option", { value: "coverage" }, "Greatest coverage"),
    el("option", { value: "name" }, "Alphabetical"));
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
    const byCoverage = order.value === "coverage";
    const peakFor = byCoverage ? peakCoverage : peak;
    const cellRank = byCoverage ? coverageRank : columnRank;
    const rank = (label) => anchor === null ? peakFor.get(label) : cellRank(hits.get(anchor).get(label));
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
          el("a", { href: contigHrefFor(row.contig_id, by) }, "open"))),
      ...shown.map((label) => pin(label, identityCell(label, hits.get(row.contig_id), otherLabel))))));
  };
  const filter = el("input", { class: "filter", type: "search", placeholder: "Filter contigs",
    "aria-label": "Filter contigs", oninput: (event) => { query = event.target.value.toLowerCase(); render(); } });
  render();
  return [el("div", { class: "controls" }, filter, el("label", {}, "Columns ", order)),
    el("div", { class: "scroll matrix" }, el("table", {}, el("thead", {}, headRow), body))];
}

// `jobHref(job)` resolves where a job's name should link to; null renders plain text instead of
// a link. Default matches the authenticated SPA's own job pages; the public SPA passes its own
// (a job only links out when that genome has also been made individually visible).
const defaultJobHref = (job) => `#/job/${job.id}`;

function compositionBars(jobs, jobHref = defaultJobHref) {
  const totals = new Map();
  for (const job of jobs) for (const row of job.composition) {
    if (row.percent_composition > 0) totals.set(row.genospecies, (totals.get(row.genospecies) ?? 0) + row.percent_composition);
  }
  if (!totals.size) return [];
  // Colour follows the species for the whole batch; anything past the seventh folds into Other --
  // but every name that landed there stays one click away, not just lost in a grey "Other" swatch.
  const named = [...totals.keys()].sort((a, b) => totals.get(b) - totals.get(a) || a.localeCompare(b)).slice(0, 7);
  const others = [...totals.keys()].filter((name) => !named.includes(name))
    .sort((a, b) => totals.get(b) - totals.get(a));
  const slot = (name) => named.includes(name) ? `var(--c${named.indexOf(name)})` : "var(--muted)";
  const swatch = (label, colour) => {
    const mark = el("span", { class: "swatch" });
    mark.style.background = colour;
    return el("span", {}, mark, label);
  };
  const otherDetail = el("div", { class: "muted other-detail", hidden: "" },
    others.map((name) => `${name} ${number(totals.get(name), 1)}%`).join(", "));
  const otherToggle = others.length ? el("button", { class: "swatch-toggle", "aria-expanded": "false",
    onclick: (event) => {
      const expanded = event.currentTarget.getAttribute("aria-expanded") === "true";
      event.currentTarget.setAttribute("aria-expanded", String(!expanded));
      otherDetail.hidden = expanded;
    } }, swatch(`Other (${others.length})`, "var(--muted)")) : null;
  const legend = el("div", {},
    el("div", { class: "legend" }, ...named.map((name) => swatch(name, slot(name))), otherToggle),
    others.length ? otherDetail : null);
  const nameOrLink = (job) => {
    const href = jobHref(job);
    return href ? el("a", { href }, job.name) : job.name;
  };
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
    return el("div", { class: "track composition" }, el("div", { class: "track-label", title: job.name }, nameOrLink(job)),
      el("div", {}, stack, el("div", { class: "muted" }, text || "No classified contigs")));
  });
  return [el("h2", {}, "Genospecies composition"), legend, ...rows];
}

function repliconTable(jobs, jobHref = defaultJobHref) {
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
  const nameOrLink = (job) => {
    const href = jobHref(job);
    return href ? el("a", { href }, job.name) : job.name;
  };
  return [el("h2", {}, "Called replicons (kb)"),
    el("p", { class: "muted" }, "Bases called to each replicon, shaded within each column. A multi-locus contig counts towards each of its names."),
    el("div", { class: "scroll matrix" }, el("table", {},
      el("thead", {}, el("tr", {}, el("th", { scope: "col" }, "Assembly"),
        ...names.map((name) => el("th", { class: "num" }, name)))),
      el("tbody", {}, ...jobs.map((job) => el("tr", {},
        el("th", { scope: "row", title: job.name }, nameOrLink(job)),
        ...names.map((name) => cell(job, name)))))))];
}

function tabs(panels, initial = 0) {
  const bodies = panels.map(([, nodes]) => el("div", {}, ...nodes));
  const buttons = panels.map(([label], index) => el("button", { class: "tab", onclick: () => show(index) }, label));
  const show = (index) => {
    bodies.forEach((body, i) => { body.hidden = i !== index; });
    buttons.forEach((button, i) => button.setAttribute("aria-pressed", i === index));
  };
  show(initial < panels.length ? initial : 0);
  return [el("div", { class: "tabs" }, ...buttons), ...bodies];
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

function overlapTrack(parts, length) {
  // One bar, every part's intervals drawn on it in its own colour, translucent so overlaps show.
  const bar = el("div", { class: "bar fusion-bar" });
  parts.forEach(({ hit }, index) => {
    if (!hit) return;
    for (const [start, end] of hit.query_intervals) {
      const span = el("span", {});
      span.style.left = `${(start - 1) / length * 100}%`;
      span.style.width = `${Math.max((end - start + 1) / length * 100, 0.2)}%`;
      span.style.background = `var(--c${index % 7})`;
      bar.append(span);
    }
  });
  return bar;
}

function fusionLegend(parts) {
  return el("div", { class: "legend" }, ...parts.map(({ label }, index) => {
    const mark = el("span", { class: "swatch" });
    mark.style.background = `var(--c${index % 7})`;
    return el("span", {}, mark, label);
  }));
}

function fusionView(tiling, plasmids, length) {
  // `tiling` is the caller's own greedy decomposition (tile_contig): the real parts it picked to
  // explain the contig, each guaranteed to have its own unfused hit -- not every fusion-named BLAST
  // hit the contig happens to score against, most of which are noise with no real standalone part.
  if (tiling.length < 2) return [el("p", { class: "muted" }, "No fusion tiling on this contig.")];
  const parts = tiling.map(([label], index) => ({ label, index, hit: plasmids.find((cell) => cell.label === label) }));
  const sideBySide = parts.map(({ label, hit }) => hit
    ? hitTrack(label, hit, length, { detail: [
        `${number(hit.overall_percent_identity, 2)}% identity over ${number(hit.query_coverage_percent, 1)}% of the contig (${number(hit.query_covered_length)} bp)`,
      ] })
    : el("div", { class: "track" }, el("div", { class: "track-label" }, label),
        el("div", { class: "muted" }, "no hit data"), el("div")));
  return [
    el("p", { class: "muted" }, `Tiling (pick order): ${tiling.map(([label, bp]) => `${label} (+${number(bp)} bp)`).join(", then ")}.`),
    el("h3", {}, "Side by side"),
    ...sideBySide,
    el("h3", {}, "Combined: overlapping on one line"),
    el("div", { class: "track" }, el("div", { class: "track-label" }, "tiling"), overlapTrack(parts, length), el("div")),
    fusionLegend(parts),
  ];
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
