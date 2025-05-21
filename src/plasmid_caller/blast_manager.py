from __future__ import annotations
import subprocess
from pathlib import Path
from importlib import resources
import sys

class BlastManager:
    """Manages BLAST binary access and execution"""

    def __init__(self) -> None:
        self.module_dir = Path(__file__).parent
        # locate manage_blast.sh as part of the package.
        with resources.as_file(
            resources.files("plasmid_caller.scripts").joinpath("manage_blast.sh")
        ) as script_path:
            self.manage_script = script_path.resolve()

        self.vendor_dir = self.manage_script.parent.parent / "vendor" / "blast"
        self.vendor_dir.mkdir(parents=True, exist_ok=True)
        self.binaries: list[str] = ["blastn", "blastp", "blastx", "tblastn", "tblastx", "makeblastdb", "blastdbcmd"]
        self._blast_path: Path | None = None
        self._blast_source: str | None = None
        self._versions: str | None = None
        # todo: convert to dict: zip - binaries | versions
        # todo: add logging

    @property
    def blast_path(self) -> Path:
        """Path object to the directory that holds BLAST binaries."""
        if self._blast_path is None:
            self._initialise_blast()
        return self._blast_path

    @property
    def blast_source(self) -> str | None:
        """get the source of the blast (system or local). This needs to use the check from the bash script."""
        if self._blast_source is None:
            _ = self.blast_path
        return self._blast_source

    @property
    def versions(self) -> str:
        """return a string where each line is a binary and its version"""
        if self._versions is None:
            _ = self._get_versions()
        return self._versions

    def _initialise_blast(self) -> None:
        """Run the shell helper and record its answer."""
        # capture via pipe and merge stderr to stdout for live display.
        proc = subprocess.Popen(
            [str(self.manage_script), "path"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        )
        last_line = ""
        assert proc.stdout is not None

        for line in proc.stdout:
            print(line, end="")
            last_line = line
        proc.wait()

        if proc.returncode != 0:
            raise RuntimeError(f"manage_blast.sh failed with exit-code {proc.returncode}")

        if not last_line.strip():
            raise RuntimeError("manage_blast.sh produced no path on stdout")

        installation_path = Path(last_line.strip())
        if not (installation_path.is_dir() and (installation_path / "blastn").exists()):
            raise RuntimeError(f"Invalid BLAST path reported: {installation_path}")

        self._blast_path = installation_path

        # Source info if available
        src_file = self.vendor_dir / "BLAST_SOURCE"
        if src_file.is_file():
            self._blast_source = src_file.read_text().strip()

    def get_binary_path(self, binary_name) -> Path:
        """get full path to a specified BLAST binary"""
        return self.blast_path / binary_name

    def run_blast_command(self, binary, *args, **kwargs):
        """run any blast command using the managed binaries."""
        print(f"Executing blast command using: {binary}")
        binary_path = self.get_binary_path(binary)
        if not binary_path.exists():
            raise FileNotFoundError(binary_path)

        cmd = [str(binary_path), *map(str, args)]
        if kwargs:
            for key, value in kwargs.items():
                if key.startswith("_"):
                    continue
                if isinstance(value, bool):
                    if value:
                        cmd.append(f"-{key}")
                else:
                    cmd.extend([f"-{key}", str(value)])
        try:
            print(f"Executing command: '{' '.join(cmd)}'")
            return subprocess.run(cmd, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"BLAST command failed ({' '.join(cmd)})\n"
                f"stdout:\n{e.stdout}\nstderr:\n{e.stderr}"
            ) from e

    def _get_versions(self) -> None:
        """set self.version to a string where each line is composed of the blast binary path and its version"""
        lines = []
        for binary in self.binaries:
            path = self.get_binary_path(binary)
            proc = subprocess.run( [str(path), "-version"], check=True, text=True, capture_output=True)
            lines.append(f"{binary}\t{proc.stdout.splitlines()[0].split(': ')[1]}")
        self._versions = "\n".join(lines)