from typing import Dict, List
import pandas as pd

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
# Toutes les autres substitutions (hors identiques et N) seront des transversions. [web:12]

def classify_mutation(ref_nt: str, alt_nt: str) -> str:
    if ref_nt == alt_nt:
        return "same"
    if (ref_nt, alt_nt) in TRANSITIONS:
        return "transition"
    return "transversion"

def find_snps_from_fasta(
    ref_seq: str,
    sequences: Dict[str, str],
) -> pd.DataFrame:
    """
    Compare chaque séquence à la référence position par position.
    Retourne un DataFrame avec colonnes:
    sequence_id, position, ref_nt, alt_nt, mutation_type
    position: index 1-based pour être cohérent avec les usages VCF. [web:9]
    """
    snp_records: List[dict] = []
    ref_len = len(ref_seq)

    for seq_id, seq in sequences.items():
        min_len = min(ref_len, len(seq))
        for i in range(min_len):
            ref_nt = ref_seq[i]
            alt_nt = seq[i]
            if ref_nt == "N" or alt_nt == "N":
                continue
            if ref_nt != alt_nt:
                mut_type = classify_mutation(ref_nt, alt_nt)
                snp_records.append(
                    {
                        "sequence_id": seq_id,
                        "position": i + 1,
                        "ref_nt": ref_nt,
                        "alt_nt": alt_nt,
                        "mutation_type": mut_type,
                    }
                )

    return pd.DataFrame(snp_records)

def find_snps_from_vcf(df_vcf: pd.DataFrame) -> pd.DataFrame:
    """
    Transforme un DataFrame VCF en DataFrame standard SNP.
    Ici on ne fait pas d'appel de variants, on suppose le VCF déjà prêt. [web:7][web:19]
    Colonnes de sortie:
    sequence_id (CHROM), position, ref_nt, alt_nt, mutation_type
    """
    records = []
    for _, row in df_vcf.iterrows():
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        ref = str(row["REF"])
        alt_field = str(row["ALT"])
        # ALT peut contenir plusieurs allèles séparés par des virgules
        alt_alleles = alt_field.split(",")

        for alt in alt_alleles:
            if len(ref) == 1 and len(alt) == 1:
                mut_type = classify_mutation(ref, alt)
                records.append(
                    {
                        "sequence_id": chrom,
                        "position": pos,
                        "ref_nt": ref,
                        "alt_nt": alt,
                        "mutation_type": mut_type,
                    }
                )
            else:
                # Pour simplifier, on ignore les INDELs ici
                continue

    return pd.DataFrame(records)
