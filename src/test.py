import pandas as pd, numpy as np
from sklearn.metrics import average_precision_score
from canonicalize import canon_pep, canon_hla

test = pd.read_csv("data/clean/itsndb_test.csv")
tp = pd.read_csv("data/clean/bigmhc_pred.csv"); tp["key"] = tp.pep.map(canon_pep)+"|"+tp.mhc.map(canon_hla)
cal = pd.read_csv("data/clean/bigmhc_calibrated.csv")
m = test.merge(tp[["key","BigMHC_IM"]], on="key").merge(cal, on="key")
non = m[m.hla != "HLA-A*02:01"].dropna(subset=["cal_z"]).reset_index(drop=True)
y, raw, calz = non.label.values, non.BigMHC_IM.values, non.cal_z.values

rng = np.random.default_rng(0); idx = np.arange(len(y)); d = []
for _ in range(5000):
    b = rng.choice(idx, len(idx), replace=True)
    if 0 < y[b].sum() < len(b):
        d.append(average_precision_score(y[b], calz[b]) - average_precision_score(y[b], raw[b]))
d = np.array(d)
print(f"non-A2 paired AUPRC delta (cal-raw): median {np.median(d):.3f}, "
      f"95% CI [{np.percentile(d,2.5):.3f}, {np.percentile(d,97.5):.3f}], P(delta>0)={(d>0).mean():.3f}")
