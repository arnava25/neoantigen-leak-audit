import pandas as pd
import mhcgnomes

AA = set("ACDEFGHIKLMNPQRSTVWY")

def canon_hla(raw):
    """Any HLA format -> 'HLA-A*02:01'; None if unparseable or not clean class-I 4-digit."""
    if pd.isna(raw):
        return None
    try:
        s = mhcgnomes.parse(str(raw).strip()).to_string()
        return s if s.startswith("HLA-") and "*" in s else None
    except Exception:
        return None

def canon_pep(raw):
    """Uppercase, strip whitespace, reject non-standard residues."""
    if pd.isna(raw):
        return None
    s = "".join(str(raw).strip().upper().split())
    return s if s and set(s) <= AA else None

def canonicalize(df, pep_col, hla_col, label_col, label_map, source,
                 keep_lengths=(8, 9, 10, 11)):
    out = pd.DataFrame({
        "peptide": df[pep_col].map(canon_pep),
        "hla":     df[hla_col].map(canon_hla),
        "label":   df[label_col].map(label_map),
        "source":  source,
    })
    out["length"] = out["peptide"].str.len()
    out["key"]    = out["peptide"] + "|" + out["hla"]   # canonical key for dedup/leakage
    n_bad = out[["peptide", "hla", "label"]].isna().any(axis=1).sum()
    out = out.dropna(subset=["peptide", "hla", "label"])
    out = out[out.length.isin(keep_lengths)].reset_index(drop=True)
    return out, int(n_bad)



if __name__ == "__main__":
    import os
    os.makedirs("data/clean", exist_ok=True)

    # --- test set: ITSNdb ---
    its = pd.read_csv("ITSNdb/data/ITSNdb.csv")
    its_c, its_bad = canonicalize(its, "Neoantigen", "HLA", "NeoType",
                                  {"Positive": 1, "Negative": 0}, "ITSNdb")
    its_c.to_csv("data/clean/itsndb_test.csv", index=False)
    print(f"ITSNdb test: {len(its_c)} rows ({its_bad} dropped), "
          f"labels {its_c.label.value_counts().to_dict()}, {its_c.hla.nunique()} alleles")

    # --- training negatives: NEPdb (HLA-I, clean 'peptide' rows only) ---
    nep = pd.read_csv("NEPdb_VND.csv")
    nep = nep[(nep.locus.isin(["A", "B", "C"])) & (nep.antigen_type == "peptide")].copy()
    nep_c, _ = canonicalize(nep, "mut_peptide", "alleleA", "response",
                            {"P": 1, "N": 0}, "NEPdb")
    before = len(nep_c)
    nep_c = nep_c[~nep_c.key.isin(set(its_c.key))].reset_index(drop=True)
    nep_c.to_csv("data/clean/nepdb_train.csv", index=False)
    print(f"NEPdb train: {len(nep_c)} rows after removing {before - len(nep_c)} "
          f"that leak into the test set, labels {nep_c.label.value_counts().to_dict()}")

    # --- primary training source: Müller harmonized neo-peptides ---
    mul = pd.read_csv("Neopep_data_org.txt", sep="\t", low_memory=False)
    print("\nMüller raw rows:", len(mul))
    print("response_type values:", mul.response_type.value_counts(dropna=False).to_dict())

    # scope = SNV only; keep ONLY screened rows: CD8 -> positive, negative -> negative.
    # 'not_tested' (the vast majority) maps to NaN via the label_map and is dropped.
    mul = mul[mul.mutation_type == "SNV"].copy()
    mul_c, _ = canonicalize(mul, "mutant_seq", "mutant_best_alleles", "response_type",
                            {"CD8": 1, "negative": 0}, "Muller")
    print("Müller SNV / screened / 8-11mer / parseable HLA:", len(mul_c),
          "labels:", mul_c.label.value_counts().to_dict())

    # mandatory leakage removal vs the test set
    before, before_pos = len(mul_c), int((mul_c.label == 1).sum())
    mul_c = mul_c[~mul_c.key.isin(set(its_c.key))].reset_index(drop=True)
    after_pos = int((mul_c.label == 1).sum())
    mul_c.to_csv("data/clean/muller_train.csv", index=False)
    print(f"Müller train after dedup: {len(mul_c)} rows "
          f"(removed {before - len(mul_c)} leaking into test; positives {before_pos} -> {after_pos}), "
          f"labels {mul_c.label.value_counts().to_dict()}")