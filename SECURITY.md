# Security Policy

GenomeX is an analysis pipeline: it reads FASTA files, shells out to
bioinformatics tools, and writes reports. It has no network listener, no
authentication, and no credential handling. The realistic risk surface is
command construction and file handling.

## Supported versions

| Version | Supported |
|---|---|
| 0.1.x (`main`) | yes |

Pre-1.0 — only the latest `main` receives fixes.

## Reporting a vulnerability

Use **GitHub private vulnerability reporting**: the *Security* tab →
*Report a vulnerability*. That opens a private channel with the maintainer.
Please do not open a public issue for anything exploitable.

Include the input that triggers it, the command, and what you observed.
A crafted FASTA file that reproduces the problem is ideal.

Expect an acknowledgement within **7 days** and an assessment within **30**.
This is a single-maintainer research project, not a funded product; those are
honest targets rather than a guarantee of a same-day patch.

## What counts

- Argument injection through file names or user-supplied parameters reaching
  `subprocess` in `genomex/runtime.py`.
- Path traversal through an input path or an output directory.
- Resource exhaustion from a malformed FASTA that is disproportionate to input
  size (an unbounded allocation, not merely a large genome).
- Anything that causes GenomeX to write outside `--outdir`.

## What does not

- Wrong biological calls. Those are correctness bugs — open a normal issue,
  and see `CONTRIBUTING.md`, where they are the most wanted contribution.
- Vulnerabilities in Prodigal, HMMER, MMseqs2, fastANI or Mash. Report those
  upstream; tell us if a workaround belongs here.
- Running GenomeX on untrusted input as root. Do not do that.
