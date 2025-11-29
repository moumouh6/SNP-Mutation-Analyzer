import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_snps_per_sequence(df_counts: pd.DataFrame):
    """
    Histogramme / barplot: nombre de SNPs par séquence.
    Retourne la figure matplotlib.
    """
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.barplot(
        data=df_counts,
        x="sequence_id",
        y="num_snps",
        ax=ax,
        color="#4C72B0",
    )
    ax.set_xlabel("Séquence")
    ax.set_ylabel("Nombre de SNPs")
    ax.set_title("Mutations par séquence")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    return fig

def plot_mutation_type_pie(df_types: pd.DataFrame):
    """
    Pie chart: proportion de types de mutations.
    """
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.pie(
        df_types["count"],
        labels=df_types["mutation_type"],
        autopct="%1.1f%%",
        startangle=90,
    )
    ax.set_title("Répartition des types de mutations")
    ax.axis("equal")
    return fig

def plot_snp_heatmap(df_snps: pd.DataFrame, ref_len: int = None):
    """
    Heatmap: fréquences des SNPs par position et par séquence.
    On construit une matrice (sequence_id x position) avec count de SNPs.
    Pour des gros génomes, il faudra limiter la taille ou binner les positions. [web:9]
    """
    if df_snps.empty:
        return None

    pivot = (
        df_snps.groupby(["sequence_id", "position"])
        .size()
        .reset_index(name="count")
        .pivot(index="sequence_id", columns="position", values="count")
        .fillna(0)
    )

    if ref_len is not None and ref_len < 2000:
        # Optionnel: s'assurer que toutes les positions existent, sinon seaborn s'en sort déjà bien.
        pass

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.heatmap(pivot, cmap="viridis", ax=ax)
    ax.set_xlabel("Position")
    ax.set_ylabel("Séquence")
    ax.set_title("Heatmap des SNPs (séquence x position)")
    plt.tight_layout()
    return fig
