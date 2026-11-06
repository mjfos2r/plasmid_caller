"use strict";
// Authenticated SPA: upload, batch/job history, visibility control. Shared renderers live in common.js.
const view = document.getElementById("view");
let poll = null;

const statusCell = (job) => el("td", { class: `status-${job.status}` }, job.status);
const ACCEPTED_NAME = /\.(fa|fna|fasta)(\.gz)?$/i;

// A dropped *folder* never shows up in dataTransfer.files -- the browser only exposes it
// through DataTransferItem.webkitGetAsEntry(), which hands back a tree of File/Directory
// entries that has to be walked by hand (and paged: a directory reader can return entries in
// batches, so it must be called repeatedly until it comes back empty). Safari has no entries
// API at all, so that case falls back to the flat file list it does give us.
function readEntries(reader) {
  return new Promise((resolve, reject) => reader.readEntries(resolve, reject));
}

async function collectEntry(entry, out) {
  if (entry.isFile) {
    const file = await new Promise((resolve, reject) => entry.file(resolve, reject));
    out.push(file);
  } else if (entry.isDirectory) {
    const reader = entry.createReader();
    let batch;
    do {
      batch = await readEntries(reader);
      for (const child of batch) await collectEntry(child, out);
    } while (batch.length);
  }
}

async function filesFromDrop(dataTransfer) {
  const entries = [...dataTransfer.items].map((item) => item.webkitGetAsEntry?.());
  if (!entries.some(Boolean)) return [...dataTransfer.files];
  const out = [];
  for (const entry of entries) if (entry) await collectEntry(entry, out);
  return out;
}

function dropZone(onFiles) {
  const input = el("input", { type: "file", multiple: "", accept: ".fa,.fna,.fasta,.gz",
    "aria-label": "Choose assembly files" });
  const folderInput = el("input", { type: "file", multiple: "", webkitdirectory: "",
    "aria-label": "Choose a folder of assemblies" });
  const zone = el("div", { class: "drop" },
    el("label", {}, input,
      el("strong", {}, "Drop assemblies or folders here"),
      el("div", { class: "muted" }, "or click to choose .fa, .fna or .fasta files, optionally gzipped")),
    el("label", { class: "drop-folder muted" }, folderInput, "or choose a folder"));
  input.addEventListener("change", () => onFiles(input.files));
  folderInput.addEventListener("change", () => onFiles(folderInput.files));
  zone.addEventListener("dragover", (event) => { event.preventDefault(); zone.classList.add("over"); });
  zone.addEventListener("dragleave", () => zone.classList.remove("over"));
  zone.addEventListener("drop", async (event) => {
    event.preventDefault();
    zone.classList.remove("over");
    onFiles(await filesFromDrop(event.dataTransfer));
  });
  return zone;
}

function uploadToggle(onFiles) {
  const button = el("button", { class: "upload-toggle", onclick: () => button.replaceWith(dropZone(onFiles)) },
    "Upload assemblies");
  return button;
}

async function upload(fileList, message) {
  const files = [...fileList].filter((file) => ACCEPTED_NAME.test(file.name));
  const skipped = fileList.length - files.length;
  if (!files.length) {
    message.className = skipped ? "error" : "muted";
    message.textContent = skipped ? "No .fa, .fna or .fasta files found." : "";
    return;
  }
  const form = new FormData();
  for (const file of files) form.append("files", file);
  message.className = "muted";
  message.textContent = `Uploading ${files.length} file(s)${skipped ? ` (${skipped} skipped)` : ""}…`;
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
  view.replaceChildren(uploadToggle((files) => upload(files, message)), message);
  const batches = await api("/api/batches");
  if (!batches.length) return;
  const rows = batches.map((batch) => {
    const done = batch.jobs.filter((job) => job.status === "done" || job.status === "found").length;
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
  const [batch, overview] = await Promise.all([api(`/api/batches/${id}`), api(`/api/batches/${id}/overview`)]);
  const viewable = (status) => status === "done" || status === "found";
  const rows = batch.jobs.map((job) => el("tr", {},
    el("td", {}, viewable(job.status) ? el("a", { href: `#/job/${job.id}` }, job.name) : job.name),
    statusCell(job),
    el("td", { class: "num" }, number(job.contigs)),
    el("td", { class: "num" }, number(job.bases)),
    el("td", {}, job.error || "")));
  view.replaceChildren(
    backButton(),
    el("h1", {}, "Batch"),
    el("p", { class: "muted" }, `Submitted ${when(batch.created_at)} · plasmid_caller ${batch.version}`),
    visibilityToggle(`/api/batches/${batch.id}/visibility`, batch.visible,
      "Visible on the public page (charts cover only the genomes you've made visible)"),
    el("div", { class: "scroll" }, el("table", {},
      el("thead", {}, el("tr", {}, el("th", {}, "Assembly"), el("th", {}, "Status"),
        el("th", { class: "num" }, "Contigs"), el("th", { class: "num" }, "Bases"), el("th", {}, "Note"))),
      el("tbody", {}, ...rows))),
    ...compositionBars(overview), ...repliconTable(overview));
  if (batch.jobs.some((job) => job.status === "queued" || job.status === "running")) {
    poll = setTimeout(() => batchPage(id).catch((error) => showError(view, error)), 2000);
  }
}

function bundleButtons(job) {
  return el("p", {},
    "Download: ",
    el("a", { href: `/api/jobs/${job.id}/download/bundle.tar.gz` }, ".tar.gz"),
    " | ",
    el("a", { href: `/api/jobs/${job.id}/download/bundle.zip` }, ".zip"),
    el("div", { class: "muted" }, `sha256 ${job.genome_hash}`));
}

function visibilityToggle(path, visible, label) {
  const checkbox = el("input", { type: "checkbox", checked: visible ? "" : null,
    onchange: (event) => api(path, {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ visible: event.target.checked }),
    }).catch((error) => { event.target.checked = !event.target.checked; alert(error.message); }) });
  return el("label", { class: "muted" }, checkbox, ` ${label}`);
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
  const contigHrefFor = (contigId, focus) => contigHref("#/job", id, contigId, focus);
  view.replaceChildren(
    backButton(),
    el("h1", {}, job.name),
    bundleButtons(job),
    job.status === "found" ? el("p", { class: "notice" }, "This genome was already called -- reusing the stored result.") : null,
    visibilityToggle(`/api/jobs/${job.id}/visibility`, job.visible, "Visible on the public page"),
    el("p", { class: "muted" }, `${number(job.contigs)} contigs · ${number(job.bases)} bp · ${summary || "no classified contigs"}`),
    ...tabs([
      ["Contigs", [
        el("h2", {}, flagged ? `Contigs (${flagged} flagged: PF32/WP disagree or ambiguous genospecies)` : "Contigs"),
        ...contigTable(contigs, contigHrefFor)]],
      ["Plasmid identity", identityTable(contigs, identity.plasmid, "Genospecies", "plasmid", contigHrefFor)],
      ["Genospecies identity", identityTable(contigs, identity.genospecies, "Plasmid", "genospecies", contigHrefFor)],
    ]),
    ...(job.provenance ? [el("p", { class: "muted" }, `plasmid_caller ${job.provenance.caller_version} · BLAST ${job.provenance.blast_version} · ` +
      `database ${job.provenance.database_sha256.slice(0, 12)} · run ${when(job.finished_at)}`)] : []),
    el("h2", {}, "Downloads"),
    el("div", { class: "links" }, ...downloads.map((name) =>
      el("a", { href: `/api/jobs/${id}/download/${name}` }, name))));
}

async function contigPage(jobId, contigId, focus) {
  const { job, contig, rule, plasmids, genospecies, pf32, loci, tiling } =
    await api(`/api/jobs/${jobId}/contig?contig_id=${encodeURIComponent(contigId)}`);
  const length = contig.contig_len;
  const initialTab = focus === "genospecies" ? 1 : 0;
  const species = genospecies.filter((cell) => cell.passes).sort((a, b) => b.score - a.score);
  const margin = species.length > 1 ? `${species[0].label} leads ${species[1].label} on coverage × identity: ` +
    `${number(species[0].score / 100, 1)} against ${number(species[1].score / 100, 1)}.` : null;
  const pf32Tracks = pf32.map((hit) => hitTrack(`${hit.plasmid_name} · ${hit.strain ?? "?"}`, hit, length, { detail: [
    `${number(hit.overall_percent_identity, 2)}% identity over ${number(hit.query_covered_length)} bp`,
    `Reference: ${hit.strain ?? "unknown strain"} ${hit.plasmid_id ?? ""}`] }));
  const plasmidHref = (cell) => el("a", { href: labelHref("#/job", jobId, "plasmid", cell.label) }, cell.label);
  const genospeciesHref = (cell) => el("a", { href: labelHref("#/job", jobId, "genospecies", cell.label) }, cell.label);
  view.replaceChildren(
    backButton(),
    el("h1", {}, contig.contig_id),
    el("p", { class: "muted" }, `${number(length)} bp · ${contig.genospecies ?? "no genospecies"}`),
    el("p", {}, "Call: ", el("strong", {}, contig.filtered_call ?? "none"), rule ? ` · ${rule}` : ""),
    copyFastaButton(`/api/jobs/${jobId}/contig/${encodeURIComponent(contigId)}/fasta`),
    fusionPanel(contig, loci, tiling),
    el("div", { class: "axis muted" }, el("span", {}, "1"), el("span", {}, `${number(length)} bp`)),
    el("h2", {}, "PF32 hits"),
    ...(pf32Tracks.length ? pf32Tracks : [el("p", { class: "muted" }, "No PF32 hits on this contig.")]),
    ...tabs([
      ["WP by plasmid", wpTracks(plasmids, length, plasmidHref)],
      ["WP by genospecies", [margin ? el("p", { class: "muted" }, margin) : null,
        ...wpTracks(genospecies, length, genospeciesHref)].filter(Boolean)],
      ...(tiling.length >= 2 ? [["Fusion view", fusionView(tiling, plasmids, length)]] : []),
    ], initialTab));
}

async function labelContigsPage(jobId, by, label) {
  const [{ job, contigs }, identity] = await Promise.all([api(`/api/jobs/${jobId}/results`), api(`/api/jobs/${jobId}/identity`)]);
  const contigLengths = new Map(contigs.map((row) => [row.contig_id, row.contig_len]));
  const cells = (by === "plasmid" ? identity.plasmid : identity.genospecies).filter((cell) => cell.label === label);
  view.replaceChildren(
    backButton(),
    el("h1", {}, label),
    ...labelContigsTable(cells, by, contigLengths, (contigId) => contigHref("#/job", jobId, contigId)));
}

async function referenceDbPage() {
  const db = await api("/api/reference-db");
  const wp = db.wp, pf32 = db.pf32;
  const wpSpecies = new Set(wp.map((row) => row.genospecies)).size;
  const wpStrains = new Set(wp.map((row) => row.strain)).size;
  const pf32Strains = new Set(pf32.map((row) => row.strain)).size;
  view.replaceChildren(
    backButton(),
    el("h1", {}, "Reference database"),
    el("p", { class: "muted" },
      `${number(wp.length)} WP references across ${number(wpStrains)} strains and ${number(wpSpecies)} genospecies · ` +
      `${number(pf32.length)} PF32 references across ${number(pf32Strains)} strains.`),
    ...tabs([
      ["WP (whole-plasmid)", sortableTable(wp, [
        { key: "genospecies", label: "Genospecies" }, { key: "strain", label: "Strain" },
        { key: "name", label: "Replicon" },
        { key: "accession", label: "Accession", render: (row) =>
          el("a", { href: ncbiNuccoreUrl(row.accession), target: "_blank", rel: "noopener" }, row.accession) },
        { key: "length", label: "Length (bp)", num: true },
        { key: "job_id", label: "Called", render: (row) => row.job_id
          ? el("a", { href: `#/job/${row.job_id}` }, "View")
          : el("span", { class: "muted" }, "not run") },
      ], "genospecies")],
      ["PF32 (ParA)", sortableTable(pf32, [
        { key: "strain", label: "Strain" }, { key: "name", label: "Replicon" },
        { key: "accession", label: "Accession" },
      ], "name")],
    ]));
}

function route() {
  clearTimeout(poll);
  const [, page, id, sub, extra, focus] = location.hash.split("/");
  const load = page === "batch" ? batchPage(id) :
    page === "job" && sub === "contig" ? contigPage(id, decodeURIComponent(extra), focus) :
    page === "job" && (sub === "plasmid" || sub === "genospecies") ? labelContigsPage(id, sub, decodeURIComponent(extra)) :
    page === "job" ? jobPage(id) :
    page === "reference-db" ? referenceDbPage() : home();
  load.catch((error) => showError(view, error));
}

window.addEventListener("hashchange", route);
route();
