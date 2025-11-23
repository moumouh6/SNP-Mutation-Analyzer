import streamlit as st

st.title("SNP & Mutation Analyzer")

uploaded_ref = st.file_uploader("Séquence de référence (FASTA)", type="fasta")
uploaded_sequences = st.file_uploader("Séquences à analyser (FASTA)", type="fasta", accept_multiple_files=True)
