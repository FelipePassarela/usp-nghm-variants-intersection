import argparse
import re
from pathlib import Path

import pandas as pd
import pysam


def main():
    parser = argparse.ArgumentParser(
        description="Find intersection of variants from VCF coming from USP and selected features from CSV."
    )
    parser.add_argument("vcf_path", type=Path, help="Path to the VCF file.")
    parser.add_argument("csv_path", type=Path, help="Path to the CSV file.")

    args = parser.parse_args()
    vcf_path = args.vcf_path
    csv_path = args.csv_path

    usp_df = usp_variants(vcf_path)

    estimator_df = estimator_variants(csv_path, n_features=20)
    intersected_df = pd.merge(
        estimator_df,
        usp_df,
        how="inner",
        left_on=["CHROM", "POS"],
        right_on=["CHROM", "POS"],
    )
    intersected_df = intersected_df[["CHROM", "POS"]].reset_index(drop=True)

    print("Intersected Variants:")
    print(intersected_df)


def usp_variants(vcf_path: Path) -> pd.DataFrame:
    vcf = pysam.VariantFile(vcf_path)
    vcf_df = pd.DataFrame(columns=["CHROM", "POS"])

    count = 0
    for rec in vcf.fetch():
        vcf_df.loc[count] = [rec.chrom, rec.pos]
        count += 1
    return vcf_df


def estimator_variants(csv_path: Path, n_features: int = 20) -> pd.DataFrame:
    """
    Read the CSV, find the row with n_features == n_features,
    and extract the variants listed in `selected_features`.
    Returns DataFrame with columns: CHROM, POS (POS as int).
    """
    df = pd.read_csv(csv_path, dtype=str)

    if "n_features" not in df.columns:
        raise KeyError("Column 'n_features' not found in CSV.")
    mask = pd.to_numeric(df["n_features"], errors="coerce") == int(n_features)

    if not mask.any():
        return pd.DataFrame(columns=["CHROM", "POS"])

    row = df.loc[mask].iloc[0]
    sel = row.get("selected_features", "") or ""

    # regex to capture tokens like 'chr14_23967207'
    matches = re.findall(r"(chr[^'\",\]\s]+?_\d+)", str(sel))

    records = []
    for m in matches:
        try:
            chrom, pos = m.split("_", 1)
            records.append({"CHROM": chrom, "POS": int(pos)})
        except Exception:
            continue

    return pd.DataFrame(records)


if __name__ == "__main__":
    main()
