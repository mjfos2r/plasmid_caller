#scripts/build_hooks.py

import os
import shutil
import stat
import subprocess
from pathlib import Path
from typing import Any, Dict

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    """
    custon build hook to handle blast installation and database setup.
    """
    def initialize(self, version: str, build_data: Dict[str, Any]) -> None:
        """Initialize the build hook."""
        print("Initialize hook called")
        print(f"Version: {version}")
        print(f"Build data keys: {list(build_data.keys() if build_data else [])}")
        return None

    def finalize(self, version: str, build_data: Dict[str, Any], artifact_path: str) -> None:
        """run after building the package to ensure script permissions and DB setup."""
        directory = Path(artifact_path).parent if artifact_path else Path.cwd()
        print(f"Using directory: {directory}")
        script_path = directory / "scripts" / "manage_blast.sh"
        if script_path.exists():
            # make it executable in case it isn't
            current_mode = os.stat(script_path).st_mode
            os.chmod(script_path, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        # make sure vendors dir exists.
        vendor_dir = Path(directory) / "vendor" / "blast"
        vendor_dir.mkdir(parents=True, exist_ok=True)

        # check for blast and install if we don't have it.
        try:
            subprocess.run(
                [str(script_path), "path"],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            print("BLAST is available")
        except subprocess.CalledProcessError:
            print("Failed to find BLAST. Installation required.")
            subprocess.run([str(script_path), "install"], check=False)
        except Exception as e:
            print(f"Error checking BLAST: {e}")
