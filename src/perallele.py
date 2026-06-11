import pandas as pd
from sklearn.metrics import roc_auc_score
from canonicalize import canon_pep, canon_hla

t = pd.read_csv("data/clean/itsndb_test.csv")
p = pd.read_csv("data/clean/bigmhc_pred.csv")
p["key"] = p["pep"].map(canon_pep) + "|" + p["mhc"].map(canon_hla)
m = t.merge(p[["key", "BigMHC_IM"]], on="key").rename(columns={"BigMHC_IM": "score"})

nonA2 = m[m.hla != "HLA-A*02:01"]
rows = []
for hla, g in nonA2.groupby("hla"):
    y = g.label.values
    if len(g) >= 6 and 0 < y.sum() < len(y):
        rows.append((hla, len(g), int(y.sum()), round(roc_auc_score(y, g.score), 3)))
rep = pd.DataFrame(rows, columns=["hla", "n", "pos", "auroc"]).sort_values("n", ascending=False)
print("per-allele AUROC (non-A*02:01, n>=6):")
print(rep.to_string(index=False) if len(rep) else "  (no non-A*02:01 allele has n>=6)")

rare = nonA2.groupby("hla").filter(lambda g: len(g) < 6)
y = rare.label.values
if 0 < y.sum() < len(y):
    print(f"\nrare alleles pooled (n<6 each): n={len(rare)} pos={int(y.sum())} "
          f"AUROC={roc_auc_score(y, rare.score):.3f}")
