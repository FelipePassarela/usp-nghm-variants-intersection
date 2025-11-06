import argparse
from pathlib import Path

import pandas as pd
import pysam


def main():
    parser = argparse.ArgumentParser(
        description="Show genotypes frequencies from a VCF file."
    )
    parser.add_argument("vcf_path", type=Path, help="Path to the VCF file.")
    parser.add_argument(
        "--output", "-o", type=Path, help="Output CSV file for genotype frequencies."
    )

    args = parser.parse_args()
    vcf_path = args.vcf_path
    output_path = args.output

    df_genotypes = genotypes_matrix(vcf_path)
    df_frequencies = genotype_frequencies(df_genotypes)

    if output_path:
        df_frequencies.to_csv(output_path)
        print(f"✅ Genotype frequencies saved to {output_path}")
    else:
        print(df_frequencies)


def genotypes_matrix(vcf_path: Path) -> pd.DataFrame:
    """Read VCF and return a DataFrame with samples as rows and variants as columns.

    Create a DataFrame containing the genotypes for each sample (rows) and each variant (columns).
    Genotypes are represented as strings (i.e., "0/0", "0/1", "1/1", "./.").

    Args:
        vcf_path (Path): Path to the VCF file.

    Returns:
        pd.DataFrame: DataFrame with samples as rows and variants as columns.
    """
    vcf = pysam.VariantFile(str(vcf_path))
    samples = list(vcf.header.samples)
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
                phased = getattr(sample, "phased", False)
                sep = "|" if phased else "/"
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
    freqs = freqs.reindex(columns=["0/0", "0/1", "1/1", "./."], level=0)
    freqs.columns.name = "genotype"
    return freqs


if __name__ == "__main__":
    main()
