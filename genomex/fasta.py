"""FASTA parsing and assembly statistics. Pure Python, no external tools."""

from __future__ import annotations

import gzip
import hashlib
import io
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

_NUC = set("ACGTacgt")


def _open_maybe_gzip(path: Path) -> io.TextIOBase:
    if str(path).endswith((".gz", ".bgz")):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (header, sequence) pairs. Header excludes the leading '>'."""
    path = Path(path)
    header: str | None = None
    chunks: list[str] = []
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_fasta(path: str | Path, records: list[tuple[str, str]], width: int = 70) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")
    return path


@dataclass
class Contig:
    name: str
    seq: str
    description: str = ""

    @property
    def length(self) -> int:
        return len(self.seq)

    @property
    def gc(self) -> float:
        """GC fraction over unambiguous bases; 0.0 if the contig is all-N."""
        s = self.seq.upper()
        at = s.count("A") + s.count("T")
        gc = s.count("G") + s.count("C")
        total = at + gc
        return gc / total if total else 0.0

    @property
    def n_count(self) -> int:
        return len(self.seq) - sum(1 for c in self.seq if c in _NUC)


@dataclass
class Assembly:
    path: Path
    name: str
    contigs: list[Contig] = field(default_factory=list)

    @classmethod
    def load(cls, path: str | Path, name: str | None = None) -> "Assembly":
        path = Path(path)
        contigs = []
        for header, seq in iter_fasta(path):
            parts = header.split(None, 1)
            contigs.append(
                Contig(name=parts[0], seq=seq, description=parts[1] if len(parts) > 1 else "")
            )
        stem = path.name
        for suffix in (".gz", ".fna", ".fa", ".fasta", ".fas"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
        return cls(path=path, name=name or stem, contigs=contigs)

    def __len__(self) -> int:
        return len(self.contigs)

    @property
    def total_bp(self) -> int:
        return sum(c.length for c in self.contigs)

    def sha256(self) -> str:
        h = hashlib.sha256()
        with open(self.path, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
        return h.hexdigest()


def _nx(lengths: list[int], total: int, fraction: float) -> tuple[int, int]:
    """Return (Nx, Lx): contig length at which cumulative coverage passes `fraction`."""
    acc = 0
    for i, length in enumerate(sorted(lengths, reverse=True), start=1):
        acc += length
        if acc >= total * fraction:
            return length, i
    return 0, 0


def assembly_stats(asm: Assembly) -> dict:
    lengths = [c.length for c in asm.contigs]
    total = sum(lengths)
    if not lengths:
        return {"n_contigs": 0, "total_bp": 0}
    n50, l50 = _nx(lengths, total, 0.5)
    n90, l90 = _nx(lengths, total, 0.9)
    gc_bases = sum(c.seq.upper().count("G") + c.seq.upper().count("C") for c in asm.contigs)
    acgt = sum(
        sum(c.seq.upper().count(b) for b in "ACGT") for c in asm.contigs
    )
    return {
        "n_contigs": len(lengths),
        "total_bp": total,
        "longest_contig": max(lengths),
        "shortest_contig": min(lengths),
        "mean_contig": round(total / len(lengths), 1),
        "n50": n50,
        "l50": l50,
        "n90": n90,
        "l90": l90,
        "gc_percent": round(100.0 * gc_bases / acgt, 2) if acgt else 0.0,
        "n_bases": sum(c.n_count for c in asm.contigs),
        "contigs_ge_1kb": sum(1 for x in lengths if x >= 1000),
        "contigs_ge_10kb": sum(1 for x in lengths if x >= 10_000),
    }
