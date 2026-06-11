import pandas as pd, numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score
from canonicalize import canon_pep, canon_hla

tr = pd.read_csv("data/clean/im_train.csv")
tr["key"] = tr.pep.map(canon_pep) + "|" + tr.mhc.map(canon_hla)
train_keys = set(tr.key.dropna())

test = pd.read_csv("data/clean/itsndb_test.csv")
pp = pd.read_csv("data/clean/prime_pred.csv")
pp["key"] = pp.peptide.map(canon_pep) + "|" + pp.hla.map(canon_hla)
m = test.merge(pp[["key", "prime_score"]], on="key")

leaked = m.key.isin(train_keys)
print(f"leaked: {leaked.sum()}/{len(m)} (positives {int(m.loc[leaked,'label'].sum())}/{int(m.label.sum())})")
clean = m[~leaked]

def block(df, name):
    y, s = df.label.values, df.prime_score.values
    if y.sum() in (0, len(y)): print(f"[{name}] n={len(df)} degenerate"); return
    print(f"[{name}] n={len(df)} pos={int(y.sum())} prev={y.mean():.3f}  "
          f"raw PRIME AUPRC={average_precision_score(y,s):.3f} AUROC={roc_auc_score(y,s):.3f} "
          f"P@10={y[np.argsort(-s)[:10]].mean():.2f} P@20={y[np.argsort(-s)[:20]].mean():.2f}")

print("\n=== PRIME raw, leakage-removed ===")
block(clean, "OVERALL")
block(clean[clean.hla == "HLA-A*02:01"], "A*02:01")
block(clean[clean.hla != "HLA-A*02:01"], "non-A*02:01")
