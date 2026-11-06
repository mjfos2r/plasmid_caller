"use strict";
// Public, unauthenticated, read-only. No upload or visibility-toggle code reachable from here --
// not hidden behind a flag, just absent, so this file can never mutate anything.
const view = document.getElementById("view");

async function publicHome() {
  const [genomes, batches] = await Promise.all([api("/api/public/genomes"), api("/api/public/batches")]);
  view.replaceChildren(el("h1", {}, "plasmid_caller: published results"));
  if (batches.length) {
    const batchRows = batches.map((batch) => {
      const n = batch.jobs.filter((job) => job.status === "done" || job.status === "found").length;
      return el("tr", {},
        el("td", {}, el("a", { href: `#/batch/${batch.id}` }, when(batch.created_at))),
        el("td", { class: "num" }, number(n)));
    });
    view.append(el("h2", {}, "Batches"),
      el("div", { class: "scroll" }, el("table", {},
        el("thead", {}, el("tr", {}, el("th", {}, "Submitted"), el("th", { class: "num" }, "Assemblies"))),
        el("tbody", {}, ...batchRows))));
  }
  if (!genomes.length) {
    view.append(el("p", { class: "muted" }, "No individual genomes are published yet."));
    return;
  }
  view.append(el("h2", {}, "Genomes"));
  const rows = genomes.map((genome) => el("tr", {},
    el("td", {}, el("a", { href: `#/genome/${genome.hash}` }, genome.name)),
    el("td", {}, genome.species ?? ""),
    el("td", {}, genome.strain ?? ""),
    el("td", {}, genome.accession_id
      ? el("a", { href: ncbiAssemblyUrl(genome.accession_id), target: "_blank", rel: "noopener" }, genome.accession_id)
      : ""),
    el("td", { class: "num" }, number(genome.contigs)),
    el("td", { class: "num" }, number(genome.bases)),
    el("td", {}, when(genome.created_at))));
  view.append(el("div", { class: "scroll" }, el("table", {},
    el("thead", {}, el("tr", {}, el("th", {}, "Assembly"), el("th", {}, "Species"), el("th", {}, "Strain"),
      el("th", {}, "Accession"), el("th", { class: "num" }, "Contigs"),
      el("th", { class: "num" }, "Bases"), el("th", {}, "Updated"))),
    el("tbody", {}, ...rows))));
}

async function genomePage(hash) {
  const [{ genome, contigs, composition, downloads }, identity] = await Promise.all([
    api(`/api/public/genomes/${hash}`), api(`/api/public/genomes/${hash}/identity`)]);
  for (const row of contigs) {
    const signals = [row.multiple_loci_pf32 && "PF32", identity.candidates[row.contig_id] && "WP"].filter(Boolean);
    row.fusion_signal = signals.length ? signals.join(" + ") : null;
  }
  const classified = composition.filter((row) => row.genospecies !== "unclassified");
  const summary = classified.map((row) => `${row.genospecies} ${number(row.percent_composition, 1)}%`).join(" · ");
  const flagged = contigs.filter(isFlagged).length;
  const contigHrefFor = (contigId, focus) => contigHref("#/genome", hash, contigId, focus);
  view.replaceChildren(
    backButton(),
    el("h1", {}, genome.name),
    el("p", {},
      "Download: ",
      el("a", { href: `/api/public/genomes/${hash}/download/bundle.tar.gz` }, ".tar.gz"),
      " | ",
      el("a", { href: `/api/public/genomes/${hash}/download/bundle.zip` }, ".zip"),
      el("div", { class: "muted" }, `sha256 ${hash}`)),
    el("p", { class: "muted" }, `${number(genome.contigs)} contigs · ${number(genome.bases)} bp · ${summary || "no classified contigs"}`),
    ...tabs([
      ["Contigs", [
        el("h2", {}, flagged ? `Contigs (${flagged} flagged: PF32/WP disagree or ambiguous genospecies)` : "Contigs"),
        ...contigTable(contigs, contigHrefFor)]],
      ["Plasmid identity", identityTable(contigs, identity.plasmid, "Genospecies", "plasmid", contigHrefFor)],
      ["Genospecies identity", identityTable(contigs, identity.genospecies, "Plasmid", "genospecies", contigHrefFor)],
    ]),
    ...(genome.provenance ? [el("p", { class: "muted" },
      `plasmid_caller ${genome.provenance.caller_version} · BLAST ${genome.provenance.blast_version} · ` +
      `database ${genome.provenance.database_sha256.slice(0, 12)} · updated ${when(genome.created_at)}`)] : []),
    el("h2", {}, "Downloads"),
    el("div", { class: "links" }, ...downloads.map((name) =>
      el("a", { href: `/api/public/genomes/${hash}/download/${name}` }, name))));
}

async function genomeContigPage(hash, contigId, focus) {
  const { genome, contig, rule, plasmids, genospecies, pf32, loci, tiling } =
    await api(`/api/public/genomes/${hash}/contig?contig_id=${encodeURIComponent(contigId)}`);
  const length = contig.contig_len;
  const initialTab = focus === "genospecies" ? 1 : 0;
  const species = genospecies.filter((cell) => cell.passes).sort((a, b) => b.score - a.score);
  const margin = species.length > 1 ? `${species[0].label} leads ${species[1].label} on coverage × identity: ` +
    `${number(species[0].score / 100, 1)} against ${number(species[1].score / 100, 1)}.` : null;
  const pf32Tracks = pf32.map((hit) => hitTrack(`${hit.plasmid_name} · ${hit.strain ?? "?"}`, hit, length, { detail: [
    `${number(hit.overall_percent_identity, 2)}% identity over ${number(hit.query_covered_length)} bp`,
    `Reference: ${hit.strain ?? "unknown strain"} ${hit.plasmid_id ?? ""}`] }));
  const plasmidHref = (cell) => el("a", { href: labelHref("#/genome", hash, "plasmid", cell.label) }, cell.label);
  const genospeciesHref = (cell) => el("a", { href: labelHref("#/genome", hash, "genospecies", cell.label) }, cell.label);
  view.replaceChildren(
    backButton(),
    el("h1", {}, contig.contig_id),
    el("p", { class: "muted" }, `${number(length)} bp · ${contig.genospecies ?? "no genospecies"}`),
    el("p", {}, "Call: ", el("strong", {}, contig.filtered_call ?? "none"), rule ? ` · ${rule}` : ""),
    copyFastaButton(`/api/public/genomes/${hash}/contig/${encodeURIComponent(contigId)}/fasta`),
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

async function labelContigsPage(hash, by, label) {
  const [{ genome, contigs }, identity] = await Promise.all([api(`/api/public/genomes/${hash}`), api(`/api/public/genomes/${hash}/identity`)]);
  const contigLengths = new Map(contigs.map((row) => [row.contig_id, row.contig_len]));
  const cells = (by === "plasmid" ? identity.plasmid : identity.genospecies).filter((cell) => cell.label === label);
  view.replaceChildren(
    backButton(),
    el("h1", {}, label),
    ...labelContigsTable(cells, by, contigLengths, (contigId) => contigHref("#/genome", hash, contigId)));
}

async function publicBatchPage(id) {
  const [batch, overview] = await Promise.all([api(`/api/public/batches/${id}`), api(`/api/public/batches/${id}/overview`)]);
  const jobHref = (job) => job.visible ? `#/genome/${job.genome_hash}` : null;
  const rows = overview.map((job) => {
    const href = jobHref(job);
    return el("tr", {},
      el("td", {}, href ? el("a", { href }, job.name) : job.name),
      el("td", {}, job.species ?? ""),
      el("td", {}, job.strain ?? ""),
      el("td", {}, job.accession_id
        ? el("a", { href: ncbiAssemblyUrl(job.accession_id), target: "_blank", rel: "noopener" }, job.accession_id)
        : ""),
      el("td", { class: "num" }, number(job.contigs)),
      el("td", { class: "num" }, number(job.bases)));
  });
  view.replaceChildren(
    backButton(),
    el("h1", {}, "Batch"),
    el("p", { class: "muted" }, `Submitted ${when(batch.created_at)} · plasmid_caller ${batch.version} · ${number(overview.length)} assemblies`),
    el("div", { class: "scroll" }, el("table", {},
      el("thead", {}, el("tr", {}, el("th", {}, "Assembly"), el("th", {}, "Species"), el("th", {}, "Strain"),
        el("th", {}, "Accession"), el("th", { class: "num" }, "Contigs"), el("th", { class: "num" }, "Bases"))),
      el("tbody", {}, ...rows))),
    ...compositionBars(overview, jobHref), ...repliconTable(overview, jobHref));
}

function route() {
  const [, page, id, sub, extra, focus] = location.hash.split("/");
  const load = page === "genome" && sub === "contig" ? genomeContigPage(id, decodeURIComponent(extra), focus) :
    page === "genome" && (sub === "plasmid" || sub === "genospecies") ? labelContigsPage(id, sub, decodeURIComponent(extra)) :
    page === "genome" ? genomePage(id) :
    page === "batch" ? publicBatchPage(id) : publicHome();
  load.catch((error) => showError(view, error));
}

window.addEventListener("hashchange", route);
route();
