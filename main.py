"""Show genotypes frequencies from a VCF file, separated by cohorts."""

import argparse
import re
from pathlib import Path

import pandas as pd
import pysam


def main():
    parser = argparse.ArgumentParser(
        description="Show genotypes frequencies from a VCF file."
    )
    parser.add_argument("vcf_path", type=Path, help="Path to the VCF file.")
    parser.add_argument(
        "cohorts",
        nargs="+",
        type=Path,
        help="Paths to cohort CSV files (space separated).",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=Path,
        help="Output directory for frequency CSV files.",
    )

    args = parser.parse_args()
    vcf_path = args.vcf_path
    cohorts_path = args.cohorts
    output_dir: Path = args.output_dir

    cohort_ids = get_cohort_ids(cohorts_path)

    for cohort_name in cohort_ids.columns:
        print(f"⚙️ Processing cohort: {cohort_name}")

        ids_to_keep = cohort_ids[cohort_name].dropna().tolist()
        df_genotypes = genotypes_matrix(vcf_path, ids_to_keep=ids_to_keep)
        df_frequencies = genotype_frequencies(df_genotypes)

        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            cohort_output_path = output_dir / f"{cohort_name}.csv"
            df_frequencies.to_csv(cohort_output_path)
            print(f"✅ Genotype frequencies saved to {cohort_output_path}")
        else:
            print(df_frequencies)


def get_cohort_ids(cohorts_path: list[Path]) -> pd.DataFrame:
    """Get cohort IDs from CSV files.

    Args:
        cohorts_path (list[Path]): List of paths to cohort CSV files.

    Returns:
        pd.DataFrame: DataFrame with cohort IDs as columns.

    Raises:
        ValueError: If no matching columns are found in a cohort.
    """

    cohorts = [
        pd.read_csv(path, sep=None, engine="python", header=None)
        for path in cohorts_path
    ]
    max_len = max(len(cohort) for cohort in cohorts)

    for i, cohort in enumerate(cohorts):
        # The files do not have a consistent header, so we search for the ID column
        id_pattern = re.compile(r"^C\d+-ExC\d+-xgenV\d+$")
        col_matches = [
            col
            for col in cohort.columns
            if cohort[col].astype(str).str.match(id_pattern).any()
        ]
        id_col = next(iter(col_matches), None)
        if id_col is None:
            raise ValueError(
                f"No matching {id_pattern.pattern} columns found in cohort."
            )

        # Only ID column is required
        cohort = cohort[[id_col]].dropna().reset_index(drop=True)
        cohorts[i] = cohort

        # Some cohorts have fewer IDs; pad with NaNs
        if len(cohort) < max_len:
            diff = max_len - len(cohort)
            padding = pd.DataFrame({id_col: [pd.NA] * diff})
            cohorts[i] = pd.concat([cohort, padding], ignore_index=True)

    cohort_ids = pd.concat(cohorts, axis=1)
    cohort_ids.columns = [path.stem for path in cohorts_path]
    return cohort_ids


def genotypes_matrix(
    vcf_path: Path, ids_to_keep: list[str] | None = None
) -> pd.DataFrame:
    """Read VCF and return a DataFrame with samples as rows and variants as columns.

    Create a DataFrame containing the genotypes for each sample (rows) and each variant (columns).
    Genotypes are represented as strings (i.e., "0/0", "0/1", "1/1", "./.").

    Args:
        vcf_path (Path): Path to the VCF file.
        ids_to_keep (list[str]): List of sample IDs to include in the DataFrame. If None, include all samples.

    Returns:
        pd.DataFrame: DataFrame with samples as rows and variants as columns.
    """
    vcf = pysam.VariantFile(str(vcf_path))
    samples = list(vcf.header.samples)
    if ids_to_keep is not None:
        ids_set = set(ids_to_keep)  # For O(1) lookups
        samples = [s for s in samples if s in ids_set]

    data = {s: [] for s in samples}
    variant_ids = []

    for rec in vcf.fetch():
        vid = f"{rec.chrom}_{rec.pos}"
        variant_ids.append(vid)
        for s in samples:
            sample = rec.samples[s]
            gt = sample.get("GT")

            if gt is None:
                gstr = "./."
            else:
                sep = "/"
                # None alleles as '.'
                gstr = sep.join("." if a is None else str(a) for a in gt)
            data[s].append(gstr)

    df = pd.DataFrame(data, index=variant_ids).T
    df.index.name = "sample"
    df.columns.name = "variant"
    return df


def genotype_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    """Genotype frequencies for each variant.

    Given a DataFrame with samples as rows and variants as columns,
    return a DataFrame with variants as rows and genotype frequencies as columns.

    Args:
        df (pd.DataFrame): DataFrame with samples as rows and variants as columns.

    Returns:
        pd.DataFrame: DataFrame with variants as rows and genotype frequencies as columns.
    """
    counts = df.apply(pd.Series.value_counts).fillna(0).astype(int).T
    totals = counts.sum(axis=1)
    freqs = counts.div(totals, axis=0)

    # Sort columns according to typical genotype order: hom ref > het > hom alt > missing
    freqs = freqs.reindex(columns=["0/0", "0/1", "1/1", "./."], level=0)
    freqs.columns.name = "genotype"

    return freqs


if __name__ == "__main__":
    main()
