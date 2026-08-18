"""External-tool execution with provenance capture.

GenomeX assumes it runs on Linux (WSL) with the bioinformatics tools on PATH --
the Windows side is a cockpit that shells in, nothing scientific executes there.
Every external invocation is recorded: argv, exit code, wall time, tool version.
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path

INSTALL_HINT = (
    "micromamba create -y -n gx -c conda-forge -c bioconda "
    "prodigal hmmer mmseqs2 mash fastani"
)

VERSION_CMD = {
    "prodigal": (["prodigal", "-v"], None),
    "hmmsearch": (["hmmsearch", "-h"], None),
    "mmseqs": (["mmseqs", "version"], None),
    "mash": (["mash", "--version"], None),
    "fastANI": (["fastANI", "--version"], None),
}


class ToolMissing(RuntimeError):
    def __init__(self, tool: str):
        super().__init__(
            f"required tool not found on PATH: {tool!r}\n"
            f"  install with:  {INSTALL_HINT}\n"
            f"  then run genomex inside that environment "
            f"(micromamba run -n gx python -m genomex ...)"
        )
        self.tool = tool


@dataclass
class Invocation:
    tool: str
    argv: list[str]
    returncode: int
    seconds: float
    stdout_tail: str = ""
    stderr_tail: str = ""


@dataclass
class Runtime:
    """Resolves tools, runs them, and accumulates a provenance record."""

    log: list[Invocation] = field(default_factory=list)
    versions: dict[str, str] = field(default_factory=dict)
    threads: int = max(1, (os.cpu_count() or 4) - 1)

    def which(self, tool: str) -> str | None:
        return shutil.which(tool)

    def require(self, *tools: str) -> None:
        missing = [t for t in tools if not self.which(t)]
        if missing:
            raise ToolMissing(missing[0])
        for t in tools:
            self.version(t)

    def version(self, tool: str) -> str:
        if tool in self.versions:
            return self.versions[tool]
        cmd, _ = VERSION_CMD.get(tool, ([tool, "--version"], None))
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            text = (p.stdout + p.stderr).strip().splitlines()
            v = next((ln.strip() for ln in text if ln.strip()), "unknown")
        except Exception as exc:  # noqa: BLE001 - version probing must never be fatal
            v = f"unknown ({exc.__class__.__name__})"
        self.versions[tool] = v[:120]
        return self.versions[tool]

    def run(
        self,
        argv: list[str],
        *,
        cwd: str | Path | None = None,
        stdout_path: str | Path | None = None,
        timeout: int = 3600,
        check: bool = True,
    ) -> Invocation:
        tool = argv[0]
        if not self.which(tool):
            raise ToolMissing(tool)
        self.version(tool)
        t0 = time.time()
        out_fh = open(stdout_path, "wb") if stdout_path else subprocess.PIPE
        try:
            p = subprocess.run(
                argv,
                cwd=str(cwd) if cwd else None,
                stdout=out_fh if stdout_path else subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=timeout,
            )
        finally:
            if stdout_path:
                out_fh.close()
        inv = Invocation(
            tool=tool,
            argv=argv,
            returncode=p.returncode,
            seconds=round(time.time() - t0, 2),
            stdout_tail=("" if stdout_path else (p.stdout or b"").decode("utf-8", "replace")[-2000:]),
            stderr_tail=(p.stderr or b"").decode("utf-8", "replace")[-2000:],
        )
        self.log.append(inv)
        if check and p.returncode != 0:
            raise RuntimeError(
                f"{tool} failed (exit {p.returncode})\n"
                f"  argv: {' '.join(argv)}\n"
                f"  stderr: {inv.stderr_tail[-1500:]}"
            )
        return inv

    def provenance(self) -> dict:
        return {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "threads": self.threads,
            "tool_versions": dict(self.versions),
            "invocations": [asdict(i) for i in self.log],
        }

    def write_provenance(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.provenance(), indent=2), encoding="utf-8")
        return path
