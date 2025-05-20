from __future__ import annotations
import subprocess
from pathlib import Path
from importlib import resources
import sys

class BlastManager:
    """Manages BLAST binary access and execution"""

    def __init__(self) -> None:
        print("Initializing BlastManager instance. Please stand by...")
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

    @property
    def blast_path(self) -> Path:
        """Path object to the directory that holds BLAST binaries."""
        if self._blast_path is None:
            print("Getting path to blast installation")
            self._initialise_blast()
        return self._blast_path

    @property
    def blast_source(self) -> None | str:
        """get the source of the blast (system or local). This needs to use the check from the bash script. TODO"""
        if self._blast_source is None:
            _ = self.blast_path
        return self._blast_source

    def _initialise_blast(self) -> None:
        """Run the shell helper and record its answer."""
        try:
            print("Initializing BLAST installation. Please stand by...")
            proc = subprocess.run(
                [str(self.manage_script), "path"],
                check=True,
                text=True,
                capture_output=True,
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"manage_blast.sh failed with exit-code {e.returncode}\n"
                f"stdout:\n{e.stdout}\nstderr:\n{e.stderr}"
            ) from e

        # Use only the LAST line – manage_blast.sh DOES echo many messages throughout execution.
        path_line = proc.stdout.splitlines()[-1].strip()
        self._blast_path = Path(path_line)

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

    def get_versions(self) -> None:
        for binary in self.binaries:
            binary_path = self.get_binary_path(binary)
            proc = subprocess.run( [str(binary_path), "-version"], check=True, text=True, capture_output=True)
            print(f"{binary_path}\t{proc.stdout.splitlines()[0].split(': ')[1]}")