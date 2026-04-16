#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import shutil
import subprocess
import sys
import time
import tomllib


def sha256_file(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


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


def capture_command(command: list[str]) -> dict:
    command = resolve_command(command)
    proc = subprocess.run(command, capture_output=True, text=True)
    return {
        "command": command,
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
    }


def nonempty(value: str | None) -> bool:
    return bool(value and value.strip() and value.strip() != "REPLACE_ME")


def build_picard_command(config: dict, sample: dict, output_bam: pathlib.Path, output_metrics: pathlib.Path) -> list[str]:
    picard = config["picard"]
    command = list(picard["command"])
    command.extend(
        [
            f"INPUT={sample['input_bam']}",
            f"OUTPUT={output_bam}",
            f"METRICS_FILE={output_metrics}",
            f"VALIDATION_STRINGENCY={picard['validation_stringency']}",
            f"COMPRESSION_LEVEL={picard['compression_level']}",
            f"ADD_PG_TAG_TO_READS={str(picard['add_pg_tag_to_reads']).lower()}",
            f"CLEAR_DT={str(picard['clear_dt']).lower()}",
            f"REMOVE_DUPLICATES={str(picard['remove_duplicates']).lower()}",
            f"REMOVE_SEQUENCING_DUPLICATES={str(picard['remove_sequencing_duplicates']).lower()}",
        ]
    )
    command.extend(picard.get("extra_args", []))
    return command


def build_markdup_command(config: dict, sample: dict, output_bam: pathlib.Path, output_metrics: pathlib.Path) -> list[str]:
    markdup = config["markdup"]
    command = list(markdup["command"])
    command.extend(
        [
            sample["input_bam"],
            "-o",
            str(output_bam),
            "-M",
            str(output_metrics),
            "-@",
            str(markdup["threads"]),
        ]
    )
    assume_sort_order = markdup.get("assume_sort_order")
    if nonempty(assume_sort_order):
        command.extend(["--assume_sort_order", assume_sort_order])
    if markdup.get("remove_duplicates"):
        command.append("--remove_duplicates")
    command.extend(markdup.get("extra_args", []))
    return command


def build_compare_command(config: dict, sample_id: str, expected_bam: pathlib.Path, actual_bam: pathlib.Path, expected_metrics: pathlib.Path, actual_metrics: pathlib.Path, output_json: pathlib.Path) -> list[str]:
    command = list(config["compare"]["command"])
    command.extend(
        [
            "--sample-id",
            sample_id,
            "--expected-bam",
            str(expected_bam),
            "--actual-bam",
            str(actual_bam),
            "--expected-metrics",
            str(expected_metrics),
            "--actual-metrics",
            str(actual_metrics),
            "--output-json",
            str(output_json),
        ]
    )
    return command


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the Hetzner Picard-vs-markdup parity batch.")
    parser.add_argument("--manifest", required=True, help="Path to the parity TOML manifest.")
    parser.add_argument(
        "--artifacts-root",
        help="Override the manifest's corpus.results_root for this run.",
    )
    parser.add_argument(
        "--sample-id",
        action="append",
        dest="sample_ids",
        help="Restrict the run to one or more sample IDs from the manifest.",
    )
    args = parser.parse_args()

    manifest_path = pathlib.Path(args.manifest).resolve()
    repo_root = pathlib.Path(__file__).resolve().parents[2]
    config = tomllib.loads(manifest_path.read_text(encoding="utf-8"))

    results_root = (
        pathlib.Path(args.artifacts_root).resolve()
        if args.artifacts_root
        else (repo_root / config["corpus"]["results_root"]).resolve()
    )
    results_root.mkdir(parents=True, exist_ok=True)

    selected_ids = set(args.sample_ids or [])
    samples = [
        sample
        for sample in config["samples"]
        if not selected_ids or sample["id"] in selected_ids
    ]
    if not samples:
        raise SystemExit("No samples selected from manifest.")

    run_summary = {
        "status": "passed",
        "manifest": str(manifest_path),
        "results_root": str(results_root),
        "environment": {
            "label": config["environment"]["label"],
            "picard_version": config["picard"]["version"],
            "picard_jar": config["picard"]["jar"],
            "picard_jar_sha256": config["picard"]["jar_sha256"],
            "java_version": capture_command([config["environment"]["java"], "-version"]),
            "samtools_version": capture_command([config["environment"]["samtools"], "--version"]),
        },
        "samples": [],
    }

    for sample in samples:
        sample_dir = results_root / sample["id"]
        sample_dir.mkdir(parents=True, exist_ok=True)
        sample_summary = {
            "id": sample["id"],
            "cohort": sample.get("cohort"),
            "input_bam": sample["input_bam"],
            "input_sha256": None,
            "input_sha256_match": None,
            "baseline_mode": None,
            "steps": [],
        }

        input_bam = pathlib.Path(sample["input_bam"])
        input_sha256 = sha256_file(input_bam)
        sample_summary["input_sha256"] = input_sha256
        if nonempty(sample.get("input_sha256")):
            sample_summary["input_sha256_match"] = input_sha256 == sample["input_sha256"]
            if not sample_summary["input_sha256_match"]:
                sample_summary["steps"].append(
                    {
                        "name": "input_checksum",
                        "status": "failed",
                        "expected": sample["input_sha256"],
                        "actual": input_sha256,
                    }
                )
                sample_summary["status"] = "failed"
                run_summary["status"] = "failed"
                run_summary["samples"].append(sample_summary)
                continue

        if nonempty(sample.get("expected_picard_bam")) and nonempty(sample.get("expected_picard_metrics")):
            expected_picard_bam = pathlib.Path(sample["expected_picard_bam"])
            expected_picard_metrics = pathlib.Path(sample["expected_picard_metrics"])
            sample_summary["baseline_mode"] = "reused_expected_outputs"
        else:
            expected_picard_bam = sample_dir / "picard.bam"
            expected_picard_metrics = sample_dir / "picard.metrics.txt"
            picard_command = build_picard_command(config, sample, expected_picard_bam, expected_picard_metrics)
            result = run_command(
                picard_command,
                repo_root,
                sample_dir / "picard.stdout.txt",
                sample_dir / "picard.stderr.txt",
            )
            result["name"] = "picard"
            sample_summary["steps"].append(result)
            sample_summary["baseline_mode"] = "generated_in_run"
            if result["exit_code"] != 0:
                sample_summary["status"] = "failed"
                run_summary["status"] = "failed"
                run_summary["samples"].append(sample_summary)
                continue

        markdup_bam = sample_dir / "markdup-wea.bam"
        markdup_metrics = sample_dir / "markdup-wea.metrics.txt"
        markdup_command = build_markdup_command(config, sample, markdup_bam, markdup_metrics)
        result = run_command(
            markdup_command,
            repo_root,
            sample_dir / "markdup.stdout.txt",
            sample_dir / "markdup.stderr.txt",
        )
        result["name"] = "markdup_wea"
        sample_summary["steps"].append(result)
        if result["exit_code"] != 0:
            sample_summary["status"] = "failed"
            run_summary["status"] = "failed"
            run_summary["samples"].append(sample_summary)
            continue

        compare_json = sample_dir / "compare.json"
        compare_command = build_compare_command(
            config,
            sample["id"],
            expected_picard_bam,
            markdup_bam,
            expected_picard_metrics,
            markdup_metrics,
            compare_json,
        )
        result = run_command(
            compare_command,
            repo_root,
            sample_dir / "compare.stdout.txt",
            sample_dir / "compare.stderr.txt",
        )
        result["name"] = "compare"
        sample_summary["steps"].append(result)
        sample_summary["compare_json"] = str(compare_json)
        sample_summary["status"] = "passed" if result["exit_code"] == 0 else "failed"
        if result["exit_code"] != 0:
            run_summary["status"] = "failed"

        write_json(sample_dir / "summary.json", sample_summary)
        run_summary["samples"].append(sample_summary)

    write_json(results_root / "run_summary.json", run_summary)
    print(json.dumps(run_summary, indent=2))
    return 0 if run_summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
