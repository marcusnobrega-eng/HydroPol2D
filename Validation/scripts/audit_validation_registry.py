#!/usr/bin/env python3
"""Audit the HydroPol2D validation registry for completeness.

The script intentionally checks campaign metadata, not scientific validity.
Scientific pass/fail remains owned by each case-specific diagnostic summary.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys


REQUIRED_COLUMNS = {
    "case_id",
    "module",
    "status",
    "case_folder",
    "expected_behavior",
    "primary_metrics",
    "acceptance_threshold",
    "report_section",
    "figure_ids",
    "table_ids",
}

REPORT_READY_STATUSES = {"validated", "report_ready"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("registry", type=Path, help="Path to Validation_Cases.csv")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat warnings as failures.",
    )
    return parser.parse_args()


def is_blank(value: str | None) -> bool:
    return value is None or value.strip() == ""


def audit_row(row: dict[str, str], workspace: Path) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    case_id = row.get("case_id", "<missing case_id>")

    folder_value = row.get("case_folder", "")
    folder = workspace / folder_value
    if is_blank(folder_value):
        errors.append(f"{case_id}: missing case_folder")
    elif not folder.exists():
        errors.append(f"{case_id}: case_folder does not exist: {folder_value}")

    for column in (
        "expected_behavior",
        "primary_metrics",
        "acceptance_threshold",
        "report_section",
        "figure_ids",
        "table_ids",
    ):
        if is_blank(row.get(column)):
            errors.append(f"{case_id}: missing {column}")

    status = row.get("status", "").strip().lower()
    if status in REPORT_READY_STATUSES:
        validation_dir = folder / "Outputs" / "Validation"
        for filename in ("Mass_Balance.csv", "Metric_Summary.csv", "Pass_Fail.csv"):
            if not (validation_dir / filename).exists():
                errors.append(
                    f"{case_id}: report-ready case missing Outputs/Validation/{filename}"
                )
    elif "diagnostic" in status:
        warnings.append(f"{case_id}: diagnostic only, do not cite as validated")
    elif "audit_required" in status:
        warnings.append(f"{case_id}: implementation audit required before validation")
    elif "pending" in status or status == "planned":
        warnings.append(f"{case_id}: pending rerun or new diagnostic outputs")

    return errors, warnings


def main() -> int:
    args = parse_args()
    registry = args.registry.resolve()
    workspace = Path.cwd().resolve()

    if not registry.exists():
        print(f"ERROR: registry not found: {registry}", file=sys.stderr)
        return 2

    with registry.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            print("ERROR: registry has no header", file=sys.stderr)
            return 2

        missing_columns = sorted(REQUIRED_COLUMNS - set(reader.fieldnames))
        if missing_columns:
            print(
                "ERROR: registry missing required columns: "
                + ", ".join(missing_columns),
                file=sys.stderr,
            )
            return 2

        rows = list(reader)

    errors: list[str] = []
    warnings: list[str] = []
    seen: set[str] = set()

    for row in rows:
        case_id = row.get("case_id", "").strip()
        if not case_id:
            errors.append("row with missing case_id")
        elif case_id in seen:
            errors.append(f"{case_id}: duplicate case_id")
        else:
            seen.add(case_id)

        row_errors, row_warnings = audit_row(row, workspace)
        errors.extend(row_errors)
        warnings.extend(row_warnings)

    print(f"Registry: {registry}")
    print(f"Cases: {len(rows)}")
    print(f"Warnings: {len(warnings)}")
    print(f"Errors: {len(errors)}")

    for warning in warnings:
        print(f"WARNING: {warning}")
    for error in errors:
        print(f"ERROR: {error}")

    if errors or (args.strict and warnings):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
