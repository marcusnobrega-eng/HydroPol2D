"""Build input_catalog.yaml from catalog_layout.py and the notes of an existing workbook.

The catalog is the single source for the input workbook: every input once, with its
tab, block, label, kind, unit, default, when it is used and its hover note.
Usage: python make_catalog.py inventory.json input_catalog.yaml
"""
import json
import sys

import yaml

from catalog_layout import TABLE_TABS, EXCLUDED, READER_FOLLOWUPS, TABS

FLAG_TAIL = "\n\nSet 0 to disable and 1 to enable unless the flag description specifies another numeric option."
NEW_NOTES = {
    "groundwater_relative_tolerance": "groundwater_relative_tolerance is the relative change of groundwater storage that triggers the next lateral groundwater update (the lateral clock). Smaller values update more often and cost more. Default 0.001.",
    "groundwater_head_tolerance_m": "groundwater_head_tolerance_m is the water-table change in metres that triggers the next lateral groundwater update (the lateral clock). Smaller values update more often and cost more. Default 0.001 m.",
    "n_outlets_data": "n_outlets_data is the number of domain-edge cells used as outlets (rounded, at least 1). Outlets grow from the lowest perimeter cell unless an outlet-cell CSV is supplied. Default 1.",
    "Prescribed_Recharge_Timeseries_File": "Prescribed_Recharge_Timeseries_File is an optional file with a recharge time series that replaces the computed recharge to the water table. It is read only when flag_prescribed_recharge = 1.",
    "riverbed_thickness_m": "riverbed_thickness_m is the riverbed thickness in metres used in the river-aquifer exchange: conductance = riverbed Ksat x riverbed area / thickness. Default 0.5 m.",
    "Subgrid_DEM_path": "Subgrid_DEM_path is an optional high-resolution DEM GeoTIFF in metres for the lookup-table sub-grid. It is read only when the lookup sub-grid is on without overbank mode.",
    "deep_aquifer_tau_days": "deep_aquifer_tau_days is the deep-aquifer response time in days. Each subbasin reservoir drains as outflow = storage / tau, and the outflow returns to the surface at the cells weighted by DeepAquiferChannelWeight_path. Read only when flag_deep_aquifer = 1.",
    "flag_prescribed_recharge": "Replaces the computed recharge to the water table with the series in Prescribed_Recharge_Timeseries_File (Input maps tab). Default 0.",
}


def clean(note):
    if not note:
        return None
    note = note.replace(FLAG_TAIL, "")
    if note.startswith("Configuration cell. "):
        note = note[len("Configuration cell. "):]
    return note.strip()


def main(inventory_path, out_path):
    inv = json.load(open(inventory_path))
    known = {**inv["settings"], **inv["flags"]}
    seen, missing_notes, tabs = set(), [], []
    for tab in TABS:
        blocks = []
        for block in tab["blocks"]:
            rows = []
            for key, label, kind, unit, default, used in block["inputs"]:
                assert key not in seen, f"{key} appears twice"
                seen.add(key)
                note = clean((known.get(key) or {}).get("note")) or NEW_NOTES.get(key)
                if not note:
                    missing_notes.append(key)
                row = {"key": key, "label": label,
                       "kind": kind if isinstance(kind, str) else "choice",
                       "unit": unit, "used_when": used, "note": note}
                if not isinstance(kind, str):
                    row["options"] = list(kind)
                row["default"] = "required" if default is None else default
                rows.append(row)
            entry = {"title": block["title"], "inputs": rows}
            for opt in ("switch", "controlled_by", "style"):
                if opt in block:
                    entry[opt] = block[opt]
            blocks.append(entry)
        tab_out = {"name": tab["name"], "purpose": tab["purpose"], "blocks": blocks}
        if tab.get("wide_value"):
            tab_out["wide_value"] = True
        if tab.get("lists"):
            tab_out["lists"] = tab["lists"]
        tabs.append(tab_out)
    reader_keys = set(inv["settings"]) | set(inv["flags"])
    unplaced = sorted(k for k in reader_keys if k not in seen and k not in EXCLUDED)
    assert not unplaced, f"not placed: {unplaced}"
    assert not missing_notes, f"no note: {missing_notes}"
    catalog = {
        "title": "HydroPol2D input workbook",
        "developer": "Marcus N. Gomes Jr., PhD",
        "excluded": EXCLUDED,
        "reader_followups": READER_FOLLOWUPS,
        "tabs": tabs,
        "tables": TABLE_TABS,
    }
    with open(out_path, "w") as fh:
        yaml.safe_dump(catalog, fh, sort_keys=False, allow_unicode=True, width=100)
    print(f"{len(seen)} inputs in {sum(len(t['blocks']) for t in tabs)} blocks on {len(tabs)} tabs -> {out_path}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
