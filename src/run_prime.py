import os, subprocess, tempfile
import pandas as pd
import mhcgnomes

ROOT = os.path.abspath(".")
PRIME = os.path.join(ROOT, "PRIME", "PRIME")
MIX   = os.path.join(ROOT, "MixMHCpred", "MixMHCpred")

t = pd.read_csv("data/clean/itsndb_test.csv")
out_rows, missing = [], []

for hla, g in t.groupby("hla"):
    compact = mhcgnomes.parse(hla).compact_string()        # HLA-A*02:01 -> A0201
    fin = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    fin.write("\n".join(g.peptide) + "\n"); fin.close()
    fout = fin.name + ".out"
    subprocess.run([PRIME, "-i", fin.name, "-a", compact, "-o", fout, "-mix", MIX],
                   capture_output=True, text=True)
    if not os.path.exists(fout):
        missing.append((hla, compact)); continue
    o = pd.read_csv(fout, sep="\t", comment="#").rename(columns={
        "Peptide": "peptide", "Score_bestAllele": "prime_score", "%Rank_bestAllele": "prime_rank"})
    o["hla"] = hla
    out_rows.append(o[["peptide", "hla", "prime_score", "prime_rank"]])

res = pd.concat(out_rows, ignore_index=True)
res.to_csv("data/clean/prime_pred.csv", index=False)
print(f"PRIME scored {len(res)}/{len(t)} peptides across {res.hla.nunique()} alleles")
if missing:
    print("alleles PRIME could not score (dropped):", missing)
