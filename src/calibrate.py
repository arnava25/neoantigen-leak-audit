import pandas as pd, numpy as np
from sklearn.metrics import roc_auc_score
from canonicalize import canon_pep, canon_hla

test = pd.read_csv("data/clean/itsndb_test.csv")
tp = pd.read_csv("data/clean/bigmhc_pred.csv")
tp["key"] = tp["pep"].map(canon_pep) + "|" + tp["mhc"].map(canon_hla)
m = test.merge(tp[["key", "BigMHC_IM"]], on="key").rename(columns={"BigMHC_IM": "raw"})

ref = pd.read_csv("data/clean/bigmhc_ref_pred.csv")
dist = {a: g["BigMHC_IM"].to_numpy() for a, g in ref.groupby("mhc") if len(g) >= 20}

def pct(r):
    d = dist.get(r.hla);  return (d < r.raw).mean() if d is not None else np.nan
def z(r):
    d = dist.get(r.hla)
    return (r.raw - d.mean()) / d.std() if (d is not None and d.std() > 0) else np.nan
m["cal_pct"] = m.apply(pct, axis=1)
m["cal_z"]   = m.apply(z, axis=1)

non = m[m.hla != "HLA-A*02:01"]
def au(df, c):
    d = df.dropna(subset=[c]); return f"{roc_auc_score(d.label, d[c]):.3f}  (n={len(d)})"
print("non-A*02:01 raw:               ", au(non, "raw"), "   [oracle was 0.697]")
print("non-A*02:01 deployable %ile:   ", au(non, "cal_pct"))
print("non-A*02:01 deployable z-score:", au(non, "cal_z"))

out = m[["key", "cal_z"]].dropna()
out.to_csv("data/clean/bigmhc_calibrated.csv", index=False)
print(f"wrote {len(out)} calibrated predictions -> data/clean/bigmhc_calibrated.csv")