import os, subprocess, tempfile
import pandas as pd, numpy as np
import mhcgnomes
from sklearn.metrics import average_precision_score, roc_auc_score
from canonicalize import canon_pep, canon_hla

ROOT = os.path.abspath(".")
PRIME = os.path.join(ROOT, "PRIME", "PRIME")
MIX   = os.path.join(ROOT, "MixMHCpred", "MixMHCpred")

# 1. score the Müller-negative reference background with PRIME, per allele
ref_in = pd.read_csv("data/clean/bigmhc_ref_input.csv")   # cols: mhc, pep, tgt
rows, missing = [], []
for mhc, g in ref_in.groupby("mhc"):
    compact = mhcgnomes.parse(mhc).compact_string()
    fin = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    fin.write("\n".join(g.pep) + "\n"); fin.close()
    fout = fin.name + ".out"
    subprocess.run([PRIME, "-i", fin.name, "-a", compact, "-o", fout, "-mix", MIX],
                   capture_output=True, text=True)
    if not os.path.exists(fout):
        missing.append((mhc, compact)); continue
    o = pd.read_csv(fout, sep="\t", comment="#").rename(
        columns={"Peptide": "pep", "Score_bestAllele": "prime_score"})
    o["mhc"] = mhc
    rows.append(o[["mhc", "pep", "prime_score"]])
ref = pd.concat(rows, ignore_index=True)
print(f"PRIME-scored reference: {len(ref)} peptides across {ref.mhc.nunique()} alleles")
if missing: print("alleles PRIME could not score:", missing)

# 2. per-allele background distribution (>=20 refs)
dist = {a: g.prime_score.to_numpy() for a, g in ref.groupby("mhc") if len(g) >= 20}

# 3. z-transform PRIME test scores against that background
test = pd.read_csv("data/clean/itsndb_test.csv")
pp = pd.read_csv("data/clean/prime_pred.csv")
pp["key"] = pp.peptide.map(canon_pep) + "|" + pp.hla.map(canon_hla)
m = test.merge(pp[["key", "prime_score"]], on="key").rename(columns={"prime_score": "raw"})
def zscore(r):
    d = dist.get(r.hla)
    return (r.raw - d.mean()) / d.std() if (d is not None and d.std() > 0) else np.nan
m["cal_z"] = m.apply(zscore, axis=1)
m[["key", "cal_z"]].dropna().to_csv("data/clean/prime_calibrated.csv", index=False)

# 4. leakage-removed raw vs calibrated
tr = pd.read_csv("data/clean/im_train.csv")
tr["key"] = tr.pep.map(canon_pep) + "|" + tr.mhc.map(canon_hla)
clean = m[~m.key.isin(set(tr.key.dropna()))]

def block(df, name):
    d = df.dropna(subset=["cal_z"]); y = d.label.values
    if y.sum() in (0, len(y)): print(f"[{name}] n={len(d)} degenerate"); return
    print(f"[{name}] n={len(d)} pos={int(y.sum())} prev={y.mean():.3f}")
    for col, lab in [("raw", "raw"), ("cal_z", "cal")]:
        s = d[col].values
        print(f"   {lab}: AUPRC={average_precision_score(y,s):.3f} AUROC={roc_auc_score(y,s):.3f} "
              f"P@10={y[np.argsort(-s)[:10]].mean():.2f} P@20={y[np.argsort(-s)[:20]].mean():.2f}")

print("\n=== PRIME leakage-removed: raw vs calibrated ===")
block(clean, "OVERALL")
block(clean[clean.hla == "HLA-A*02:01"], "A*02:01")
nonc = clean[clean.hla != "HLA-A*02:01"]
block(nonc, "non-A*02:01")

d = nonc.dropna(subset=["cal_z"]); y = d.label.values; raw = d.raw.values; calz = d.cal_z.values
rng = np.random.default_rng(0); idx = np.arange(len(y)); dl = []
for _ in range(5000):
    b = rng.choice(idx, len(idx), replace=True)
    if 0 < y[b].sum() < len(b):
        dl.append(average_precision_score(y[b], calz[b]) - average_precision_score(y[b], raw[b]))
dl = np.array(dl)
print(f"\nnon-A2 paired AUPRC delta (cal-raw), leakage-removed: median {np.median(dl):.3f}, "
      f"95% CI [{np.percentile(dl,2.5):.3f},{np.percentile(dl,97.5):.3f}], P>0={(dl>0).mean():.3f}")
