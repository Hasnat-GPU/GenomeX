"""Standalone HTML report. No external assets, no network, no build step."""

from __future__ import annotations

import html
import json
from pathlib import Path

CSS = """
:root{--bg:#fbfbf9;--panel:#fff;--ink:#1a1a18;--muted:#6b6b64;--line:#e4e3dd;
--ok:#2f7d4f;--warn:#a86a12;--bad:#b3402f;--accent:#3d5a80;--chip:#f0efe9}
:root:not([data-theme="light"]){}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){
--bg:#16171a;--panel:#1e1f23;--ink:#e9e8e4;--muted:#9b9a93;--line:#2e3037;
--ok:#57b97e;--warn:#d99b3a;--bad:#e0705c;--accent:#7fa5cc;--chip:#25272c}}
:root[data-theme="dark"]{--bg:#16171a;--panel:#1e1f23;--ink:#e9e8e4;--muted:#9b9a93;
--line:#2e3037;--ok:#57b97e;--warn:#d99b3a;--bad:#e0705c;--accent:#7fa5cc;--chip:#25272c}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
.wrap{max-width:1080px;margin:0 auto;padding:40px 24px 80px}
h1{font-size:26px;margin:0 0 4px;letter-spacing:-.01em}
h2{font-size:19px;margin:40px 0 12px;padding-bottom:6px;border-bottom:1px solid var(--line)}
h3{font-size:15px;margin:22px 0 8px}
.sub{color:var(--muted);font-size:13px;margin-bottom:28px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:18px 20px;margin:14px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:14px}
.stat{background:var(--chip);border-radius:8px;padding:10px 12px}
.stat .k{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted)}
.stat .v{font-size:18px;font-weight:600;margin-top:2px;font-variant-numeric:tabular-nums}
.pill{display:inline-block;padding:2px 9px;border-radius:999px;font-size:12px;font-weight:600}
.pill.ok{background:color-mix(in srgb,var(--ok) 18%,transparent);color:var(--ok)}
.pill.warn{background:color-mix(in srgb,var(--warn) 18%,transparent);color:var(--warn)}
.pill.bad{background:color-mix(in srgb,var(--bad) 18%,transparent);color:var(--bad)}
.bar{display:flex;height:16px;border-radius:4px;overflow:hidden;margin:8px 0 4px}
.bar span{display:block}
.legend{font-size:12px;color:var(--muted);display:flex;gap:14px;flex-wrap:wrap}
.legend i{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:5px}
table{width:100%;border-collapse:collapse;font-size:13px;margin-top:8px}
th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line)}
th{color:var(--muted);font-weight:600;font-size:11px;text-transform:uppercase;letter-spacing:.05em}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.scroll{overflow-x:auto}
code{background:var(--chip);padding:1px 5px;border-radius:4px;font-size:12px}
.note{color:var(--muted);font-size:13px;margin-top:8px}
ul.reasons{margin:8px 0 0;padding-left:18px;color:var(--muted);font-size:13px}
.limits li{margin-bottom:6px;font-size:13px;color:var(--muted)}
"""


def _esc(x) -> str:
    return html.escape(str(x))


def _pill(verdict: str) -> str:
    cls = {"clean": "ok", "possible": "warn", "likely": "bad"}.get(verdict, "warn")
    label = {"clean": "single organism", "possible": "possible contamination",
             "likely": "likely contaminated", "undetermined": "undetermined"}.get(verdict, verdict)
    return f'<span class="pill {cls}">{_esc(label)}</span>'


def _usability_pill(usable) -> str:
    """Three states, because `usable` has three.

    `None` means the contamination check did not run, and rendering it as "hold"
    would claim a check found a problem. A boolean-shaped ternary here is how the
    unassessed case used to come out green.
    """
    if usable is None:
        return '<span class="pill warn">usability undetermined</span>'
    return (
        '<span class="pill ok">usable</span>' if usable
        else '<span class="pill bad">hold</span>'
    )


def _stat(k, v) -> str:
    return f'<div class="stat"><div class="k">{_esc(k)}</div><div class="v">{_esc(v)}</div></div>'


def _marker_bar(m: dict) -> str:
    n = max(1, m["markers_total"])
    seg = [
        (m["single_copy"], "var(--ok)", "single copy"),
        (m["duplicated"], "var(--warn)", "duplicated"),
        (m["fragmented"], "var(--accent)", "fragmented"),
        (m["missing"], "var(--bad)", "missing"),
    ]
    bars = "".join(
        f'<span style="width:{100 * c / n:.2f}%;background:{col}"></span>' for c, col, _ in seg if c
    )
    legend = "".join(
        f'<span><i style="background:{col}"></i>{lab} {c}</span>' for c, col, lab in seg
    )
    return f'<div class="bar">{bars}</div><div class="legend">{legend}</div>'


def _table(headers: list[str], rows: list[list[str]], numeric: set[int] | None = None) -> str:
    numeric = numeric or set()
    head = "".join(f"<th>{_esc(h)}</th>" for h in headers)
    body = "".join(
        "<tr>" + "".join(
            f'<td class="num">{_esc(c)}</td>' if i in numeric else f"<td>{_esc(c)}</td>"
            for i, c in enumerate(r)
        ) + "</tr>"
        for r in rows
    )
    return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def _not_measured_table(nm: dict) -> str:
    """Render the explicit refusals.

    A reader who sees no contamination section assumes the run was clean. This
    is the difference between "we looked and found nothing" and "we could not
    look", and it is the whole reason proteome results are rendered separately.
    """
    rows = [
        [key.replace("_", " "), v.get("status", ""), v.get("reason", "")]
        for key, v in sorted(nm.items())
    ]
    return _table(["not reported", "status", "why"], rows)


def render_html(data: dict) -> str:
    genomes = data.get("genomes", [])
    proteomes = data.get("proteomes", [])
    pairs = data.get("pairs", [])
    pg = data.get("pangenome")
    prov = data.get("provenance", {})

    parts: list[str] = []
    parts.append(f"<title>GenomeX Report</title><style>{CSS}</style>")
    parts.append('<div class="wrap">')
    parts.append("<h1>GenomeX report</h1>")
    counted = []
    if genomes:
        counted.append(f"{len(genomes)} genome(s)")
    if proteomes:
        counted.append(f"{len(proteomes)} proteome(s)")
    if pairs or not counted:
        counted.append(f"{len(pairs)} pairwise comparison(s)")
    parts.append(
        f'<div class="sub">{" &middot; ".join(counted)} '
        f'&middot; {data.get("seconds", 0)} s &middot; genomex {_esc(data.get("genomex_version", ""))}</div>'
    )

    # ---- per genome -------------------------------------------------------
    if genomes:
        parts.append("<h2>Genome quality</h2>")
    for g in genomes:
        a, m, c, q = g["assembly"], g["markers"], g["contamination"], g["quality_call"]
        parts.append('<div class="card">')
        parts.append(
            f'<h3>{_esc(g["genome"])} &nbsp; {_pill(c["verdict"])} '
            f'{_usability_pill(q["usable_for_comparative_analysis"])}</h3>'
        )
        parts.append('<div class="grid">')
        parts.append(_stat("size", f'{a["total_bp"] / 1e6:.2f} Mb'))
        parts.append(_stat("contigs", a["n_contigs"]))
        parts.append(_stat("N50", f'{a["n50"] / 1000:.1f} kb'))
        parts.append(_stat("GC", f'{a["gc_percent"]}%'))
        parts.append(_stat("genes", g["genes"]["n_genes"]))
        parts.append(_stat("coding density", f'{100 * g["genes"]["coding_density"]:.1f}%'))
        parts.append(_stat("completeness", f'{m["completeness_percent"]}%'))
        parts.append(_stat("duplicated", f'{m["duplication_percent"]}%'))
        parts.append("</div>")
        parts.append(f'<div class="note"><code>{_esc(m["busco_style_string"])}</code></div>')
        parts.append(_marker_bar(m))
        reasons = "".join(f"<li>{_esc(r)}</li>" for r in c["reasons"])
        # "Contamination evidence" over an abstention's reason list reads as
        # evidence of cleanliness. It is the reason the check was declined.
        heading = ("Contamination evidence" if c["suspect_contigs"] is not None
                   else "Why contamination was not assessed")
        parts.append(f'<h3>{heading}</h3><ul class="reasons">{reasons}</ul>')
        b = c.get("bins", {})
        if b.get("bin0") and b.get("bin1"):
            parts.append(
                _table(
                    ["composition bin", "contigs", "bases", "% of assembly", "mean GC%"],
                    [
                        [k, b[k]["contigs"], f'{b[k]["bp"]:,}', b[k]["bp_fraction_percent"],
                         b[k]["mean_gc_percent"]]
                        for k in ("bin0", "bin1")
                    ],
                    numeric={1, 2, 3, 4},
                )
            )
        top = g.get("_suspect_rows") or []
        if top:
            parts.append("<h3>Contigs carrying the evidence</h3>")
            parts.append(
                _table(
                    ["contig", "call", "length", "GC%", "TNF dist", "displaced markers", "flags"],
                    top, numeric={2, 3, 4, 5},
                )
            )
            parts.append(
                '<div class="note"><b>contaminant_candidate</b> foreign composition with '
                'displaced core markers &middot; <b>replicon_candidate</b> distinct but '
                'marker-clean, i.e. plasmid or second chromosome as likely as contaminant '
                '&middot; <b>atypical_host_region</b> distinct, but holds this assembly&rsquo;s '
                'only copies of core markers, so it is this organism&rsquo;s own chromosome '
                '&mdash; prophage or island, not contamination &middot; '
                '<b>marker_conflict</b> shares duplicated core markers while '
                'composition looks native.</div>'
            )
        parts.append("</div>")

    # ---- per proteome -----------------------------------------------------
    # Its own section, not a genome card with holes in it. Every stat here is a
    # property of the protein file itself; nothing is inferred about a genome.
    if proteomes:
        parts.append("<h2>Proteome completeness</h2>")
        parts.append(
            '<div class="note">Scored from a supplied protein FASTA. GenomeX did not '
            'call these genes, and there is no assembly behind them &mdash; so this '
            'section reports completeness and duplication, and nothing else.</div>'
        )
    for p in proteomes:
        m, s, q = p["markers"], p["proteins"], p["quality_call"]
        parts.append('<div class="card">')
        parts.append(
            f'<h3>{_esc(p["name"])} &nbsp; '
            f'<span class="pill {"ok" if q["completeness_grade"] == "high" else "warn"}">'
            f'{_esc(q["completeness_grade"])} completeness</span> '
            f'<span class="pill warn">contamination not measured</span></h3>'
        )
        parts.append('<div class="grid">')
        parts.append(_stat("proteins", f'{s["n_proteins"]:,}'))
        parts.append(_stat("mean length", f'{s["mean_protein_aa"]} aa'))
        parts.append(_stat("median length", f'{s["median_protein_aa"]} aa'))
        parts.append(_stat("completeness", f'{m["completeness_percent"]}%'))
        parts.append(_stat("duplicated", f'{m["duplication_percent"]}%'))
        parts.append("</div>")
        parts.append(f'<div class="note"><code>{_esc(m["busco_style_string"])}</code></div>')
        parts.append(_marker_bar(m))
        reasons = "".join(f"<li>{_esc(r)}</li>" for r in q["reasons"])
        parts.append(f'<h3>Assessment</h3><ul class="reasons">{reasons}</ul>')
        parts.append("<h3>Not measured, and why</h3>")
        parts.append(_not_measured_table(p.get("not_measured", {})))
        parts.append(
            '<div class="note"><b>impossible</b> the input does not carry the '
            'information &middot; <b>unimplemented</b> the information is derivable '
            'from a proteome, but this path does not derive it.</div>'
        )
        parts.append("</div>")

    # ---- pangenome --------------------------------------------------------
    if pg:
        parts.append("<h2>Pangenome</h2>")
        parts.append('<div class="card"><div class="grid">')
        parts.append(_stat("orthogroups", pg["orthogroups_total"]))
        parts.append(_stat("core", pg["core"]))
        parts.append(_stat("accessory", pg["accessory"]))
        parts.append(_stat("strain-unique", pg["strain_unique"]))
        parts.append(_stat("core share", f'{pg["core_percent"]}%'))
        parts.append("</div>")
        cl = pg["clustering"]
        parts.append(
            f'<div class="note">{_esc(cl["tool"])} at {cl["min_seq_id"]} identity, '
            f'{cl["coverage"]} coverage, across {pg["n_genomes"]} genomes.</div></div>'
        )

    # ---- pairs ------------------------------------------------------------
    if pairs:
        parts.append("<h2>Why these genomes differ</h2>")
    for p in pairs:
        ani = p["ani"]
        ani_txt = "not resolved" if ani["ani_percent"] is None else f'{ani["ani_percent"]}%'
        same = ani["same_species_by_ani"]
        same_txt = "same species" if same else ("different species" if same is False else "unresolved")
        parts.append('<div class="card">')
        parts.append(f'<h3>{_esc(p["genome_a"])} vs {_esc(p["genome_b"])}</h3>')
        parts.append('<div class="grid">')
        parts.append(_stat("ANI", ani_txt))
        parts.append(_stat("verdict", same_txt))
        parts.append(_stat("shared orthogroups", p["shared_orthogroups"]))
        parts.append(_stat(f'unique to {p["genome_a"]}', p["unique_to_a"]["total"]))
        parts.append(_stat(f'unique to {p["genome_b"]}', p["unique_to_b"]["total"]))
        parts.append("</div>")
        for key, label in (("unique_to_a", p["genome_a"]), ("unique_to_b", p["genome_b"])):
            u = p[key]
            if not u["total"]:
                continue
            rows = [[k, v, f'{100 * v / u["total"]:.1f}%'] for k, v in
                    sorted(u["by_explanation"].items(), key=lambda kv: -kv[1])]
            parts.append(f"<h3>Genes only in {_esc(label)} &mdash; why</h3>")
            parts.append(_table(["explanation", "genes", "share"], rows, numeric={1, 2}))
            null = u.get("island_null_test") or {}
            note = (
                f'{u["islands"]} genomic island(s) holding {u["genes_in_islands"]} genes; '
                f'{u["on_suspect_contigs"]} genes sit on contigs flagged as compositionally '
                f'foreign.'
            )
            if null.get("enrichment") is not None:
                verdict = (
                    "above chance" if null["informative"]
                    else "NOT distinguishable from chance -- do not read these as acquisitions"
                )
                note += (
                    f' Permutation null expects {null["genes_in_islands_expected_by_chance"]} '
                    f'island genes here, so clustering is {null["enrichment"]}x expected: {verdict}.'
                )
            parts.append(f'<div class="note">{_esc(note)}</div>')
        parts.append(f'<div class="note">{_esc(ani["note"])}</div>')
        parts.append("</div>")

    # ---- provenance and limits -------------------------------------------
    parts.append("<h2>Provenance</h2><div class='card'>")
    versions = prov.get("tool_versions", {})
    if versions:
        parts.append(_table(["tool", "version"], [[k, v] for k, v in sorted(versions.items())]))
    parts.append(
        f'<div class="note">{_esc(prov.get("platform", ""))} &middot; '
        f'Python {_esc(prov.get("python", ""))} &middot; '
        f'{len(prov.get("invocations", []))} external tool invocations recorded in '
        f'<code>provenance.json</code>.</div></div>'
    )

    parts.append("<h2>What this report does not tell you</h2><div class='card'><ul class='limits'>")
    proteome_limits = [
        "For a supplied proteome, completeness measures the proteome, not the genome. A marker "
        "missing here may be absent from the assembly or missed by whatever called the genes -- "
        "GenomeX did not call them and cannot separate the two.",
        "Marker duplication on a proteome is not a contamination signal. Without contigs there is "
        "no way to tell a second organism's copy from a redundant gene model.",
    ] if proteomes else []
    genome_limits = [
        "Contamination detection is composition-based. It sees organisms whose k-mer signature "
        "differs from the host; it is blind to a close relative of the host, and it cannot name "
        "the contaminant without a reference database.",
        "Contigs below the minimum length are not scored -- a short contaminant contig can pass "
        "unexamined.",
        "Gene functions are not assigned. A gene called unique here is unique in sequence-cluster "
        "terms, not annotated as a known function.",
        "Genomic islands are inferred from runs of consecutive unique genes, not from integrase "
        "detection, tRNA boundaries, or codon-usage models. Treat them as candidates.",
        "fastANI reports nothing below roughly 80% identity; an unresolved ANI means distant, "
        "not identical.",
    ] if genomes else []
    for lim in [
        "Completeness is a HMMER scan against BUSCO's odb10 profiles with BUSCO's own score and "
        "length cutoffs -- close to BUSCO, but not BUSCO. It does not re-predict genes per marker "
        "region, so counts can differ by a marker or two.",
        *proteome_limits,
        *genome_limits,
    ]:
        parts.append(f"<li>{_esc(lim)}</li>")
    parts.append("</ul></div>")

    parts.append("</div>")
    return "\n".join(parts)


def write_html_report(result, path: str | Path) -> Path:
    """Render a PipelineResult (or a plain dict) to standalone HTML."""
    data = result.to_dict() if hasattr(result, "to_dict") else result

    # Attach the top suspect contigs per genome for the table in the report.
    if hasattr(result, "genomes"):
        by_name = {g.name: g for g in result.genomes}
        for entry in data["genomes"]:
            g = by_name.get(entry["genome"])
            if not g:
                continue
            interesting = [c for c in g.contamination.contigs if c.call != "core"]
            top = sorted(interesting, key=lambda c: -c.suspicion)[:12]
            entry["_suspect_rows"] = [
                [c.name, c.call, f"{c.length:,}", c.gc_percent, c.tnf_distance,
                 c.duplicated_markers, "; ".join(c.flags) or "-"]
                for c in top
            ]

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_html(data), encoding="utf-8")
    return path


def report_from_json(json_path: str | Path, html_path: str | Path) -> Path:
    data = json.loads(Path(json_path).read_text(encoding="utf-8"))
    return write_html_report(data, html_path)
