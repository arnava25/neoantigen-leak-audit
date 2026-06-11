import sys
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

A0201 = "HLA-A*02:01"

def prec_rec_at_k(y, s, k):
    k = min(k, len(y))
    top = np.argsort(-s)[:k]
    tp = int(y[top].sum())
    return tp / k, (tp / int(y.sum()) if y.sum() else np.nan)

def bootstrap_auprc(y, s, n=2000, seed=0):
    rng = np.random.default_rng(seed); idx = np.arange(len(y)); out = []
    for _ in range(n):
        b = rng.choice(idx, len(idx), replace=True)
        if 0 < y[b].sum() < len(b):
            out.append(average_precision_score(y[b], s[b]))
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))) if out else (np.nan, np.nan)

def report(df, name, ks=(10, 20, 50)):
    y = df["label"].to_numpy().astype(int); s = df["score"].to_numpy().astype(float)
    n, npos = len(y), int(y.sum())
    print(f"\n[{name}]  n={n}  pos={npos}  prevalence={npos/n:.3f}" if n else f"\n[{name}] empty")
    if n == 0 or npos == 0 or npos == n:
        print("  (degenerate stratum; skipping metrics)"); return
    ap = average_precision_score(y, s); lo, hi = bootstrap_auprc(y, s)
    print(f"  AUPRC={ap:.3f}  (95% CI {lo:.3f}-{hi:.3f})   AUROC={roc_auc_score(y, s):.3f}")
    for k in ks:
        if k <= n:
            p, r = prec_rec_at_k(y, s, k); print(f"  P@{k}={p:.3f}  R@{k}={r:.3f}")

def evaluate(truth_csv, pred_csv, score_col, pep_col=None, hla_col=None,
             higher_is_immunogenic=True):
    truth = pd.read_csv(truth_csv)           # needs columns: key, label, hla
    pred = pd.read_csv(pred_csv)
    if "key" not in pred.columns:            # build key from peptide+hla via the canonicalizer
        from canonicalize import canon_pep, canon_hla
        pred["key"] = pred[pep_col].map(canon_pep) + "|" + pred[hla_col].map(canon_hla)
    pred = pred.rename(columns={score_col: "score"})[["key", "score"]].dropna()
    if not higher_is_immunogenic:
        pred["score"] = -pred["score"]
    pred = pred.sort_values("score").drop_duplicates("key", keep="last")  # one score per key
    m = truth.merge(pred, on="key", how="left")
    print(f"coverage: {m.score.notna().sum()}/{len(m)} test peptides scored ({m.score.notna().mean():.1%})")
    m = m.dropna(subset=["score"])
    report(m, "OVERALL")
    report(m[m.hla == A0201], "A*02:01")
    report(m[m.hla != A0201], "non-A*02:01")

if __name__ == "__main__":
    a = sys.argv  # usage: python src/eval_baseline.py <pred.csv> <score_col> [pep_col hla_col] [higher_is_immunogenic]
    evaluate("data/clean/itsndb_test.csv", a[1], a[2],
             *(a[3:5] if len(a) > 4 else ()),
             higher_is_immunogenic=(a[5].lower() != "false" if len(a) > 5 else True))
