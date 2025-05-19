from __future__ import annotations
import subprocess
from pathlib import Path
from importlib import resources
import sys

class BlastManager:
    """Manages BLAST binary access and execution"""

    def __init__(self):
        self.module_dir = Path(__file__).parent
        self.project_root = self._find_project_root()
        self.vendor_dir = self.project_root / "vendor" / "blast"
        self.scripts_dir = self.project_root / "scripts"
        self.manage_script = self.scripts_dir / "manage_blast.sh"

        self._blast_path = None
        self._blast_source = None
        # locate manage_blast.sh packaged alongside the project
        with resources.as_file(
            resources.files(__package__).joinpath("..", "scripts", "manage_blast.sh")
        ) as script_path:
            self.manage_script = script_path.resolve()

        self.vendor_dir = self.manage_script.parent.parent / "vendor" / "blast"
        self.vendor_dir.mkdir(parents=True, exist_ok=True)

        self._blast_path: Path | None = None
        self._blast_source: str | None = None

    @property
    def blast_path(self) -> Path:
        """Path object to the directory that holds BLAST binaries."""
        if self._blast_path is None:
            self._initialise_blast()
        return self._blast_path

    @property
    def blast_source(self):
        """get the source of the blast (system or local)."""
        if self._blast_source is None:
            _ = self.blast_path
        return self._blast_source

    def _initialise_blast(self) -> None:
        """Run the shell helper and record its answer."""
        try:
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

        # Use only the first line – manage_blast.sh may echo more later.
        path_line = proc.stdout.splitlines()[0].strip()
        self._blast_path = Path(path_line)

        # Source info if available
        src_file = self.vendor_dir / "BLAST_SOURCE"
        if src_file.is_file():
            self._blast_source = src_file.read_text().strip()

    def get_binary_path(self, binary_name):
        """get full path to a specified BLAST binary"""
        return self.blast_path / binary_name

    def run_blast_command(self, command, *args, **kwargs):
        """run any blast command using the managed binaries."""
        binary_path = self.get_binary_path(command)
        if not binary_path.exists():
            raise FileNotFoundError(f"{binary_path} not found")
        cmd = [str(binary_path), *map(str, args)]
        if kwargs:
            for key, value in kwargs.items():
                if key.startswith("_"):
                    continue #skip private params
            if isinstance(value, bool):
                if value:
                    cmd.append(f"-{key}")
            else:
                cmd.extend([f"-{key}", str(value)])
        try:
            return subprocess.run(cmd, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"BLAST command failed ({' '.join(cmd)})\n"
                f"stdout:\n{e.stdout}\nstderr:\n{e.stderr}"
            ) from e