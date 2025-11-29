import streamlit as st
import pandas as pd
import io
import matplotlib.pyplot as plt

def _fig_to_png(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    buf.seek(0)
    return buf.getvalue()


from utils.file_loader import (
    load_fasta_sequences,
    load_reference_fasta,
    load_vcf,
)
from utils.mutation_finder import (
    find_snps_from_fasta,
    find_snps_from_vcf,
)
from utils.stats import (
    snps_per_sequence,
    mutation_type_frequencies,
    most_frequent_positions,
    filter_by_frequency,
)
from utils.visualization import (
    plot_snps_per_sequence,
    plot_mutation_type_pie,
    plot_snp_heatmap,
)

st.set_page_config(page_title="SNP & Mutation Analyzer", layout="wide")

st.title("SNP & Mutation Analyzer")

st.sidebar.header("Import des données")

data_mode = st.sidebar.radio(
    "Type de données",
    ["FASTA (séquences)", "VCF (variants)"],
)

uploaded_ref = st.sidebar.file_uploader(
    "Séquence de référence (FASTA)",
    type=["fasta", "fa"],
)

uploaded_sequences = None
uploaded_vcf = None

if data_mode == "FASTA (séquences)":
    uploaded_sequences = st.sidebar.file_uploader(
        "Séquences à comparer (FASTA)",
        type=["fasta", "fa"],
    )
else:
    uploaded_vcf = st.sidebar.file_uploader(
        "Fichier VCF",
        type=["vcf"],
    )

st.sidebar.header("Filtres de fréquence")
min_freq = st.sidebar.slider("Fréquence minimale", 0.0, 1.0, 0.0, 0.01)
max_freq = st.sidebar.slider("Fréquence maximale", 0.0, 1.0, 1.0, 0.01)

st.sidebar.header("Export")
export_csv_button = st.sidebar.button("Exporter SNPs en CSV")
export_png_hist_button = st.sidebar.button("Exporter histogramme PNG")
export_png_pie_button = st.sidebar.button("Exporter pie chart PNG")
export_png_heatmap_button = st.sidebar.button("Exporter heatmap PNG")

df_snps = pd.DataFrame()
ref_id = None
ref_seq = None

if uploaded_ref is not None:
    ref_id, ref_seq = load_reference_fasta(uploaded_ref)
    st.success(f"Référence chargée: {ref_id} (longueur: {len(ref_seq)} nt)")
else:
    st.info("Veuillez d'abord téléverser une séquence de référence.")

if ref_seq:
    if data_mode == "FASTA (séquences)" and uploaded_sequences is not None:
        sequences = load_fasta_sequences(uploaded_sequences)
        st.write(f"Nombre de séquences: {len(sequences)}")
        if len(sequences) > 0:
            avg_len = sum(len(s) for s in sequences.values()) / len(sequences)
            st.write(f"Longueur moyenne: {avg_len:.1f} nt")
        df_snps = find_snps_from_fasta(ref_seq, sequences)

    elif data_mode == "VCF (variants)" and uploaded_vcf is not None:
        df_vcf = load_vcf(uploaded_vcf)
        st.write(f"Nombre de variants dans le VCF: {len(df_vcf)}")
        df_snps = find_snps_from_vcf(df_vcf)

if not df_snps.empty:
    st.subheader("SNPs détectés")

    # Filtrage par fréquence
    df_snps_filtered = filter_by_frequency(
        df_snps, min_freq=min_freq or None, max_freq=max_freq or None
    )

    st.dataframe(df_snps_filtered.head(500))

    col1, col2 = st.columns(2)
    with col1:
        st.markdown("### Statistiques générales")
        df_counts = snps_per_sequence(df_snps_filtered)
        df_types = mutation_type_frequencies(df_snps_filtered)
        df_top_pos = most_frequent_positions(df_snps_filtered, top_n=20)

        st.write("Nombre de SNPs par séquence")
        st.dataframe(df_counts)
        st.write("Types de mutations")
        st.dataframe(df_types)
        st.write("Positions les plus fréquemment mutées")
        st.dataframe(df_top_pos)

    with col2:
        st.markdown("### Visualisations")
        fig_hist = plot_snps_per_sequence(df_counts)
        st.pyplot(fig_hist)

        fig_pie = plot_mutation_type_pie(df_types)
        st.pyplot(fig_pie)

        fig_heatmap = plot_snp_heatmap(df_snps_filtered, ref_len=len(ref_seq))
        if fig_heatmap:
            st.pyplot(fig_heatmap)

    # Export CSV
    if export_csv_button:
        csv_data = df_snps_filtered.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="Télécharger SNPs (CSV)",
            data=csv_data,
            file_name="snps_detectes.csv",
            mime="text/csv",
        )

    # Export PNG
    if export_png_hist_button and 'fig_hist' in locals():
        st.download_button(
            label="Télécharger histogramme (PNG)",
            data=_fig_to_png(fig_hist),
            file_name="hist_snps_par_sequence.png",
            mime="image/png",
        )

    if export_png_pie_button and 'fig_pie' in locals():
        st.download_button(
            label="Télécharger pie chart (PNG)",
            data=_fig_to_png(fig_pie),
            file_name="pie_types_mutations.png",
            mime="image/png",
        )

    if export_png_heatmap_button and 'fig_heatmap' in locals() and fig_heatmap:
        st.download_button(
            label="Télécharger heatmap (PNG)",
            data=_fig_to_png(fig_heatmap),
            file_name="heatmap_snps.png",
            mime="image/png",
        )
else:
    st.info("Aucun SNP à afficher pour le moment. Charge les fichiers pour commencer.")

