"""Shared utility functions for the shufflon annotation pipeline."""

import logging
import os
import subprocess
import sys

logger = logging.getLogger("shufflonfinder")


def setup_logging(verbosity: int = 1) -> None:
    """Configure logging for the pipeline.

    Args:
        verbosity: 0 = WARNING, 1 = INFO, 2 = DEBUG
    """
    level = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}.get(
        verbosity, logging.DEBUG
    )
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )


def run_cmd(cmd: list[str], description: str = "", **kwargs) -> subprocess.CompletedProcess:
    """Run a shell command, log it, and raise on failure.

    Args:
        cmd: Command and arguments as a list of strings.
        description: Human-readable label for log messages.

    Returns:
        The CompletedProcess object.

    Raises:
        subprocess.CalledProcessError: If the command exits non-zero.
    """
    label = description or " ".join(cmd[:3])
    logger.info("Running: %s", label)
    logger.debug("Full command: %s", " ".join(cmd))
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        **kwargs,
    )
    if result.returncode != 0:
        logger.error("Command failed (%s):\nstdout: %s\nstderr: %s",
                      label, result.stdout, result.stderr)
        result.check_returncode()
    return result


def ensure_dir(path: str) -> str:
    """Create directory (and parents) if it doesn't exist. Returns the path."""
    os.makedirs(path, exist_ok=True)
    return path


def check_tool(name: str) -> str:
    """Verify an external tool is on PATH and return its full path.

    Raises:
        FileNotFoundError: If the tool is not found.
    """
    from shutil import which

    path = which(name)
    if path is None:
        raise FileNotFoundError(
            f"Required tool '{name}' not found on PATH. "
            f"Install it or activate the conda environment."
        )
    logger.debug("Found %s at %s", name, path)
    return path
