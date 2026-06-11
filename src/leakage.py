import pandas as pd, numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score
from canonicalize import canon_pep, canon_hla

tr = pd.read_csv("data/clean/im_train.csv")
tr["key"] = tr["pep"].map(canon_pep) + "|" + tr["mhc"].map(canon_hla)
train_keys = set(tr["key"].dropna())
print(f"BigMHC im_train: {len(tr)} rows, {len(train_keys)} unique peptide|HLA keys")

test = pd.read_csv("data/clean/itsndb_test.csv")
tp = pd.read_csv("data/clean/bigmhc_pred.csv"); tp["key"] = tp.pep.map(canon_pep)+"|"+tp.mhc.map(canon_hla)
cal = pd.read_csv("data/clean/bigmhc_calibrated.csv")
m = test.merge(tp[["key","BigMHC_IM"]], on="key").merge(cal, on="key")

leaked = m["key"].isin(train_keys)
print(f"test peptides in BigMHC im_train: {leaked.sum()}/{len(m)} "
      f"(positives leaked: {int(m.loc[leaked,'label'].sum())}/{int(m.label.sum())})")
clean = m[~leaked]

def block(df, name):
    d = df.dropna(subset=["cal_z"]); y = d.label.values
    if y.sum() in (0, len(y)): print(f"[{name}] n={len(d)} degenerate"); return
    print(f"[{name}] n={len(d)} pos={int(y.sum())} prev={y.mean():.3f}")
    for col, lab in [("BigMHC_IM","raw"), ("cal_z","cal")]:
        s = d[col].values
        print(f"   {lab}: AUPRC={average_precision_score(y,s):.3f} AUROC={roc_auc_score(y,s):.3f} "
              f"P@10={y[np.argsort(-s)[:10]].mean():.2f} P@20={y[np.argsort(-s)[:20]].mean():.2f}")

print("\n=== leakage-removed test ===")
block(clean, "OVERALL")
block(clean[clean.hla=="HLA-A*02:01"], "A*02:01")
nonc = clean[clean.hla!="HLA-A*02:01"]
block(nonc, "non-A*02:01")

d = nonc.dropna(subset=["cal_z"]); y=d.label.values; raw=d.BigMHC_IM.values; calz=d.cal_z.values
rng=np.random.default_rng(0); idx=np.arange(len(y)); dl=[]
for _ in range(5000):
    b=rng.choice(idx,len(idx),replace=True)
    if 0<y[b].sum()<len(b): dl.append(average_precision_score(y[b],calz[b])-average_precision_score(y[b],raw[b]))
dl=np.array(dl)
print(f"\nnon-A2 paired AUPRC delta (cal-raw), leakage-removed: median {np.median(dl):.3f}, "
      f"95% CI [{np.percentile(dl,2.5):.3f},{np.percentile(dl,97.5):.3f}], P>0={(dl>0).mean():.3f}")
