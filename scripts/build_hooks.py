#scripts/build_hooks.py

import os
import stat
import subprocess
from pathlib import Path

def post_build(directory):
    """run after package build to ensure proper permissions for our scripts"""
    script_path = Path(directory) / "scripts" / "manage_blast.sh"
    # make it executable in case it isn't
    if script_path.exists():
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
        print("Failed to find BLAST. Installation required at runtime")
    except Exception as E:
        print(f"Error checking BLAST: {e}")

    return True

