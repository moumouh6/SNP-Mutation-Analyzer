import pandas as pd

def snps_per_sequence(df_snps: pd.DataFrame) -> pd.DataFrame:
    """
    Retourne un DataFrame avec nombre total de SNPs par séquence_id.
    """
    return (
        df_snps.groupby("sequence_id")
        .size()
        .reset_index(name="num_snps")
        .sort_values("num_snps", ascending=False)
    )

def mutation_type_frequencies(df_snps: pd.DataFrame) -> pd.DataFrame:
    """
    Compte transitions vs transversions.
    """
    return (
        df_snps["mutation_type"]
        .value_counts()
        .reset_index()
        .rename(columns={"index": "mutation_type", "mutation_type": "count"})
    )

def most_frequent_positions(df_snps: pd.DataFrame, top_n: int = 20) -> pd.DataFrame:
    """
    Identifie les positions les plus fréquemment mutées (toutes séquences confondues).
    """
    return (
        df_snps.groupby("position")
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        .head(top_n)
    )

def filter_by_frequency(
    df_snps: pd.DataFrame, min_freq: float = None, max_freq: float = None
) -> pd.DataFrame:
    """
    Filtre les SNPs selon une fréquence d'apparition par position.
    La fréquence est calculée comme count(position) / total_sequences.
    À connecter avec un paramètre Streamlit (slider). [web:9]
    """
    if df_snps.empty:
        return df_snps

    counts = df_snps.groupby("position").size().reset_index(name="count")
    total_sequences = df_snps["sequence_id"].nunique()
    counts["freq"] = counts["count"] / total_sequences

    df = df_snps.merge(counts[["position", "freq"]], on="position", how="left")

    if min_freq is not None:
        df = df[df["freq"] >= min_freq]
    if max_freq is not None:
        df = df[df["freq"] <= max_freq]

    return df
