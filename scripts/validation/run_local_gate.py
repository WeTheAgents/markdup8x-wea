#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import subprocess
import sys
import time


def resolve_command(command: list[str]) -> list[str]:
    if not command:
        return command
    executable = shutil.which(command[0])
    if executable is None:
        cargo_bin = pathlib.Path.home() / ".cargo" / "bin" / command[0]
        if sys.platform.startswith("win"):
            cargo_bin = cargo_bin.with_suffix(".exe")
        if cargo_bin.exists():
            executable = str(cargo_bin)
    if executable is None:
        raise FileNotFoundError(f"Could not resolve executable: {command[0]}")
    return [executable, *command[1:]]


def run_command(command: list[str], cwd: pathlib.Path, stdout_path: pathlib.Path, stderr_path: pathlib.Path) -> dict:
    started_at = time.time()
    command = resolve_command(command)
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        proc = subprocess.run(command, cwd=cwd, stdout=stdout, stderr=stderr, text=True)
    return {
        "command": command,
        "exit_code": proc.returncode,
        "duration_seconds": round(time.time() - started_at, 3),
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the local Picard-parity proof subset.")
    parser.add_argument(
        "--manifest",
        help="Optional TOML manifest for a full parity batch. When omitted, only local Rust gates run.",
    )
    parser.add_argument(
        "--artifacts-root",
        default="validation-results/local-gate",
        help="Directory for command logs and JSON summaries.",
    )
    args = parser.parse_args()

    repo_root = pathlib.Path(__file__).resolve().parents[2]
    artifacts_root = repo_root / args.artifacts_root
    artifacts_root.mkdir(parents=True, exist_ok=True)

    summary = {
        "status": "passed",
        "repo_root": str(repo_root),
        "artifacts_root": str(artifacts_root),
        "steps": [],
    }

    steps = [
        ("cargo_test", ["cargo", "test"]),
        ("cargo_clippy", ["cargo", "clippy", "--all-targets", "--all-features", "--", "-D", "warnings"]),
    ]

    for name, command in steps:
        result = run_command(
            command,
            repo_root,
            artifacts_root / f"{name}.stdout.txt",
            artifacts_root / f"{name}.stderr.txt",
        )
        result["name"] = name
        summary["steps"].append(result)
        if result["exit_code"] != 0:
            summary["status"] = "failed"
            break

    if summary["status"] == "passed" and args.manifest:
        parity_command = [
            sys.executable,
            str(pathlib.Path(__file__).with_name("run_parity_batch.py")),
            "--manifest",
            args.manifest,
            "--artifacts-root",
            str(artifacts_root / "parity"),
        ]
        result = run_command(
            parity_command,
            repo_root,
            artifacts_root / "parity.stdout.txt",
            artifacts_root / "parity.stderr.txt",
        )
        result["name"] = "parity_batch"
        summary["steps"].append(result)
        if result["exit_code"] != 0:
            summary["status"] = "failed"
    else:
        summary["steps"].append(
            {
                "name": "parity_batch",
                "status": "skipped",
                "reason": "No manifest supplied; real-data differential validation was not run.",
            }
        )

    summary_path = artifacts_root / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
