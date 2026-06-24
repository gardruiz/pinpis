from pathlib import Path
from typing import Dict, List, Optional, Union
import os
import shutil
import subprocess


class CommandError(Exception):
    pass


def run_command(
    command: List[str],
    cwd: Optional[Union[str, Path]] = None,
    env: Optional[Dict[str, str]] = None,
) -> None:
    """
    Run an external command.

    Parameters
    ----------
    command:
        Command and arguments as a list.
    cwd:
        Optional working directory.
    env:
        Optional environment variables to add or override.
    """
    printable = " ".join(str(x) for x in command)
    print(f"[pinpis] Running: {printable}")

    run_env = os.environ.copy()

    if env:
        # subprocess expects environment values as strings
        run_env.update({str(k): str(v) for k, v in env.items()})

    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        env=run_env,
    )

    if result.returncode != 0:
        raise CommandError(
            f"Command failed with exit code {result.returncode}: {printable}"
        )


def check_executable(name: str) -> bool:
    return shutil.which(name) is not None


def require_executable(name: str) -> None:
    if not check_executable(name):
        raise CommandError(f"Required executable not found in PATH: {name}")


def ensure_dir(path: Union[str, Path]) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path
