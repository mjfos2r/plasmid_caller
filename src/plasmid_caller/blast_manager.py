import subprocess
from pathlib import Path

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

    def _find_project_root(self):
        """figure out the root directory containing vendor/scripts"""
        module_parent = self.module_dir.parent
        # are we running from src in dev mode?
        if (module_parent.parent / "scripts").exists():
            return module_parent.parent

        # are we running as an installed package?
        if (module_parent / "scripts").exists():
            return module_parent

        # if we can't find it, just use the parent dir.
        return module_parent

    @property
    def blast_path(self):
        """Get path to BLAST binaries. Install if necessary"""
        if self._blast_path is None:
            # make sure we have the vendor dir.
            self.vendor_dir.mkdir(parents=True, exist_ok=True)

            # check for the management script
            if not self.manage_script.exists():
                raise FileNotFoundError(
                    f"BLAST management script not found at {self.manage_script}."
                    "Please reinstall the package :("
                )

            # run the script and get the path
            try:
                result = subprocess.run(
                    [str(self.manage_script), "path"],
                    check=True,
                    text=True,
                    capture_output=True
                )
                self._blast_path = Path(result.stdout.strip())

                # read the source to determine if it's system or local.
                source_file = self.vendor_dir / "BLAST_SOURCE"
                if source_file.exists():
                    with open(source_file, "r") as f:
                        self._blast_source = f.read().strip()
            except subprocess.CalledProcessError as e:
                raise RuntimeError(
                    f"Failed to run BLAST management script: {e}\n"
                    f"stdout: {e.stdout}\nstderr: {e.stderr}"
                )
        return self._blast_path

    @property
    def blast_source(self):
        """get the source of the blast (system or local)."""
        if self._blast_source is None:
            _ = self.blast_path
        return self._blast_source

    def get_binary_path(self, binary_name):
        """get full path to a specified BLAST binary"""
        return self.blast_path / binary_name

    def run_blast_command(self, command, *args, **kwargs):
        """run any blast command using the managed binaries."""
        binary_path = self.get_binary_path(command)
        cmd = [str(binary_path)]
        cmd.extend(str(arg) for arg in args)
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
                f"Error running {command}: {e}\nCommand: {' '.join(cmd)}\n"
                f"stdout: {e.stdout}\nstderr: {e.stderr}"
            )
        except Exception as E:
            raise RuntimeError(f"Unexpected error running {command}: {e}")


blast_manager = BlastManager()