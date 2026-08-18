"""End-to-end test against the real toolchain.

Skipped automatically when Prodigal/HMMER/MMseqs2/fastANI are not on PATH, so the
pure-logic suite still runs anywhere. When they are present this proves the
plumbing -- argv construction, output parsing, provenance capture, report
rendering -- not just the algorithms.
"""

import json
import shutil
from pathlib import Path

import pytest

from fixtures.synthetic import clean_genome, contaminated_genome
from genomex.fasta import write_fasta
from genomex.pipeline import DEFAULT_LINEAGE, run_pipeline
from genomex.report import write_html_report

REQUIRED = ["prodigal", "hmmsearch", "mmseqs", "fastANI"]
missing = [t for t in REQUIRED if not shutil.which(t)]

pytestmark = [
    pytest.mark.tools,
    pytest.mark.skipif(bool(missing), reason=f"missing tools: {missing}"),
    pytest.mark.skipif(
        not Path(DEFAULT_LINEAGE).is_dir(), reason=f"no lineage at {DEFAULT_LINEAGE}"
    ),
]


@pytest.fixture(scope="module")
def genomes(tmp_path_factory):
    d = tmp_path_factory.mktemp("genomes")
    a = clean_genome("alpha", n_contigs=8, contig_len=25_000, seed=21)
    b = contaminated_genome(
        "beta", n_host=8, n_foreign=2, contig_len=25_000, seed=21, foreign_gc=0.35
    )
    return (
        write_fasta(d / "alpha.fna", a.records),
        write_fasta(d / "beta.fna", b.records),
        b.foreign_contigs,
    )


def test_full_pipeline_runs_and_reports(genomes, tmp_path):
    alpha, beta, planted = genomes
    result = run_pipeline([alpha, beta], tmp_path / "run", log=lambda *a: None)

    assert len(result.genomes) == 2
    assert len(result.pairs) == 1

    for g in result.genomes:
        assert g.stats["total_bp"] > 0
        assert g.calls.n_genes > 0, "Prodigal returned no genes"
        assert g.markers.n_markers == 124
        assert Path(g.outdir / "contigs.tsv").exists()

    # The planted contaminant must survive the whole pipeline, not just the unit test.
    beta_result = next(g for g in result.genomes if g.name == "beta")
    assert beta_result.contamination.suspect_contig_names() >= planted

    pair = result.pairs[0]
    assert pair.shared_orthogroups >= 0
    assert pair.ani.method == "fastANI"

    # Every external call is recorded with a version and an exit code.
    prov = result.runtime.provenance()
    assert {"prodigal", "hmmsearch"} <= set(prov["tool_versions"])
    assert all(i["returncode"] == 0 for i in prov["invocations"])
    assert any(i["tool"] == "prodigal" for i in prov["invocations"])

    json_path = result.write_json(tmp_path / "run" / "genomex.json")
    data = json.loads(json_path.read_text())
    assert data["genomes"][0]["markers"]["markers_total"] == 124

    html = write_html_report(result, tmp_path / "run" / "report.html")
    text = html.read_text(encoding="utf-8")
    assert "<title>GenomeX Report</title>" in text
    assert "What this report does not tell you" in text, "limits section must ship"
    assert "alpha" in text and "beta" in text


def test_rerun_reuses_gene_calls(genomes, tmp_path):
    """Second run into the same outdir must not re-invoke Prodigal."""
    alpha, beta, _ = genomes
    out = tmp_path / "cached"
    first = run_pipeline([alpha], out, pangenome=False, log=lambda *a: None)
    second = run_pipeline([alpha], out, pangenome=False, log=lambda *a: None)

    assert any(i["tool"] == "prodigal" for i in first.runtime.provenance()["invocations"])
    assert not any(
        i["tool"] == "prodigal" for i in second.runtime.provenance()["invocations"]
    )
    assert first.genomes[0].calls.n_genes == second.genomes[0].calls.n_genes
