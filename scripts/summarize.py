# Summarize locus genotypes per sample from Tracy BCF outputs.
# Uses ab1_manifest.csv (from sanitize step) to label QC vs samples and to resolve sanitized names.
import os, re, csv, sys, yaml
import pandas as pd
from cyvcf2 import VCF

# --- Snakemake I/O and params ---
OUT         = snakemake.output["csv"]
PANEL_YAML  = snakemake.input["panel"]
MANIFEST    = snakemake.input.get("manifest") or snakemake.params.get("manifest","")
OUTDIR      = snakemake.params["outdir"]
AB1_DIR     = snakemake.params["ab1_dir"]           # not strictly needed now, but kept for parity
SEL         = snakemake.params["loci"].split(",") if snakemake.params["loci"] else []
SHEET       = snakemake.params.get("sample_sheet") or ""

# --- Panel metadata ---
panel = yaml.safe_load(open(PANEL_YAML, "r"))
lmeta = {l["locus_id"]: l for l in panel["loci"]}

def gt_to_zyg(gt: str):
    if gt in {"0/0","0|0"}: return "hom-ref"
    if gt in {"0/1","1/0","0|1","1|0"}: return "het"
    if gt in {"1/1","1|1"}: return "hom-alt"
    return "ambiguous"

# --- Read manifest to drive sample list and control/group labels ---
manifest_rows = []
manifest_map = {}  # sanitized -> dict(row)
if MANIFEST and os.path.exists(MANIFEST):
    with open(MANIFEST, "r") as f:
        r = csv.DictReader(f)
        for row in r:
            manifest_rows.append(row)
            manifest_map[row["sanitized"]] = dict(row)
else:
    # Defensive: no manifest -> write empty CSV and exit
    pd.DataFrame([]).to_csv(OUT, index=False)
    print(f"Wrote {OUT} rows=0 loci={','.join(SEL)} (manifest missing)")
    sys.exit(0)

# Build sample list (sanitized names) from manifest
samples = []
for row in manifest_rows:
    s = row["sanitized"]
    ctrl = row["control_type"]          # "negative"|"positive"|"unknown"|"none"
    group = row["group"]                # "qc"|"samples"
    samples.append((s, ctrl, group))

records = []
vcf_dir = os.path.join(OUTDIR, "vcf")

for sbase, ctrl, group in samples:
    # Initialize row
    row = {
        "sanitized": sbase,
        "qc_group": group,                 # "qc" or "samples"
        "control_type": ctrl,              # "negative"/"positive"/"unknown"/"none"
        "qc_status": "NA",
        "qc_notes":  ""
    }

    # Initialize per-locus columns
    for locus in SEL:
        row[f"{locus}_genotype"]          = "."
        row[f"{locus}_zygosity"]          = "ambiguous"
        row[f"{locus}_mutation_present"]  = "UNKNOWN"
        row[f"{locus}_status"]            = "no_call"

    # Per-locus genotype extraction
    for locus in SEL:
        meta = lmeta.get(locus)
        if not meta:
            row["qc_notes"] += f" missing metadata for {locus};"
            continue

        bcf_path = os.path.join(vcf_dir, f"{sbase}_{locus}.bcf")

        # Missing or placeholder BCF (0 bytes) -> no_alignment
        if not os.path.exists(bcf_path) or os.path.getsize(bcf_path) == 0:
            row[f"{locus}_status"] = "no_alignment"
            continue

        try:
            v = VCF(bcf_path)
        except Exception as e:
            row["qc_notes"] += f" failed to open {locus} BCF: {e};"
            row[f"{locus}_status"] = "no_alignment"
            continue

        # Prefer rsID, else contig:pos region
        target = None
        if meta.get("rsid"):
            for rec in v:
                if (rec.ID or "").lower() == meta["rsid"].lower():
                    target = rec
                    break
            # Re-open for region query if needed later
            v.close()
            v = VCF(bcf_path)

        if target is None and meta.get("contig") and meta.get("local_pos"):
            region = f"{meta['contig']}:{meta['local_pos']}-{meta['local_pos']}"
            for rec in v(region):
                target = rec
                break

        if target:
            gt_arr = target.genotypes[0] if target.genotypes else None
            gt_str = f"{gt_arr[0]}/{gt_arr[1]}" if gt_arr and None not in gt_arr[:2] else "."
            zyg   = gt_to_zyg(gt_str)
            mut   = "YES" if zyg in ("het","hom-alt") else ("NO" if zyg=="hom-ref" else "UNKNOWN")
            row[f"{locus}_genotype"]         = gt_str
            row[f"{locus}_zygosity"]         = zyg
            row[f"{locus}_mutation_present"] = mut
            row[f"{locus}_status"]           = "called"
        else:
            # Alignment likely OK but no record at target position -> assume hom-ref
            row[f"{locus}_genotype"]         = "0/0"
            row[f"{locus}_zygosity"]         = "hom-ref"
            row[f"{locus}_mutation_present"] = "NO"
            row[f"{locus}_status"]           = "no_variant_but_good_alignment"

    # --- QC evaluation ---
    if group == "qc":
        # Both NEG and POS expected: no HFE target variant present
        any_mut = any(row.get(f"{l}_mutation_present") == "YES" for l in SEL)
        row["qc_status"] = "FAIL" if any_mut else "PASS"
        if ctrl == "negative":
            row["qc_notes"] = ("NTC shows target variant signal"
                               if any_mut else "NTC OK (no target signal)")
        elif ctrl == "positive":
            row["qc_notes"] = ("POS control shows target variant (unexpected)"
                               if any_mut else "POS control OK (no HFE target)")
        else:
            row["qc_notes"] = ("QC sample shows target variant (unexpected)"
                               if any_mut else "QC sample OK (no HFE target)")

    records.append(row)

# --- Build DataFrame and sort robustly ---
df = pd.DataFrame(records)

if df.empty:
    df.to_csv(OUT, index=False)
    print(f"Wrote {OUT} rows=0 loci={','.join(SEL)}")
    sys.exit(0)

# Prefer modern schema (qc_group + sanitized); fall back gently if needed.
if {"qc_group", "sanitized"}.issubset(df.columns):
    df["qc_rank"] = df["qc_group"].map({"qc": 0, "samples": 1}).fillna(1)
    df = df.sort_values(by=["qc_rank", "sanitized"]).drop(columns=["qc_rank"])
elif {"well", "sample_id"}.issubset(df.columns):
    # Back‑compat: older filename-parsing layout
    def well_key(w):
        if isinstance(w, str) and re.match(r"^[A-H]\d{2}$", w):
            return (w[0], int(w[1:]))
        return ("Z", 999)
    df = df.sort_values(by=["well","sample_id"],
                        key=lambda c: c.map(well_key) if c.name=="well" else c)
else:
    # Last resort: sort by any available stable id
    by_cols = [c for c in ["sanitized","sample_id"] if c in df.columns] or df.columns.tolist()
    df = df.sort_values(by=by_cols)

# Write output
df.to_csv(OUT, index=False)
print(f"Wrote {OUT} rows={len(df)} loci={','.join(SEL)}")
