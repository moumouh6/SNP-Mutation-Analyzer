from typing import List, Dict, Tuple
from Bio import SeqIO  # pip install biopython
import pandas as pd

VALID_NUCS = {"A", "C", "G", "T", "N"}

def load_fasta_sequences(path: str) -> Dict[str, str]:
    """
    Lit un fichier FASTA et retourne un dict {id_sequence: sequence}.
    """
    sequences = {}
    for record in SeqIO.parse(path, "fasta"):
        seq_str = str(record.seq).upper()
        # Filtrer éventuellement les caractères non ACGTN
        seq_str = "".join([n for n in seq_str if n in VALID_NUCS])
        sequences[record.id] = seq_str
    return sequences

def load_reference_fasta(path: str) -> Tuple[str, str]:
    """
    Lit une séquence de référence FASTA (on suppose une seule entrée).
    Retourne (ref_id, ref_sequence).
    """
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) == 0:
        raise ValueError("Aucune séquence dans le fichier de référence.")
    record = records[0]
    seq_str = str(record.seq).upper()
    seq_str = "".join([n for n in seq_str if n in VALID_NUCS])
    return record.id, seq_str

def load_vcf(path: str) -> pd.DataFrame:
    """
    Lit un VCF simple en DataFrame.
    Pour des VCF complexes ou volumineux, tu pourras basculer vers cyvcf2 plus tard. [web:7][web:19]
    """
    # On ignore les lignes de header qui commencent par '#'
    df = pd.read_csv(
        path,
        comment="#",
        sep="\t",
        header=None,
        names=[
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "SAMPLES",
        ],
    )
    return df
