# NeoForge

Leakage-controlled re-evaluation of HLA class I neoantigen immunogenicity
predictors (CD8+, single-substitution neoepitopes).

**Status: exploratory / negative result — not a manuscript.** The central
finding (leading immunogenicity tools are inflated by training-data leakage on
the standard benchmark and generalize poorly off HLA-A*02:01) was independently
published in stronger form while this was in progress — see *Related work*.
Preserved for the reproducible pipeline and data-quality notes, not as a novel
claim.

## Key findings (ITSNdb, n=199)
- **Benchmark contamination.** ~48% of ITSNdb test peptides (57% of positives)
  are in BigMHC's immunogenicity training set; PRIME shows comparable overlap
  against the shared public training pool. ITSNdb is clean as a *dataset* but
  contaminated as a *benchmark* for tools trained on the same literature.
- **Leakage-free collapse off A*02:01.** With leaked peptides removed, both
  tools hold on A*02:01 (AUROC ~0.78-0.84) but fall to chance-or-below on
  non-A*02:01 (AUROC ~0.35-0.44). The cross-allele competence was largely
  memorization.
- **Cross-allele miscalibration is tool-specific.** Per-allele recalibration
  against a presented-peptide background significantly improves BigMHC's
  non-A*02:01 ranking (paired AUPRC delta +0.167) but degrades PRIME's. No
  universal post-hoc fix; raw immunogenicity scores should not be pooled across
  alleles for ranking.

## Caveats
Leakage-free non-A*02:01 has only ~35 positives (underpowered, wide CIs);
PRIME's leakage set is approximated by BigMHC's training pool (optimistic for
PRIME); single benchmark, two tools.

## Reusable
- `src/canonicalize.py` — peptide+HLA canonicalization (mhcgnomes) + the
  data-quality filters (minimal-epitope "unit trap"; positive-leakage dedup).
- `src/eval_baseline.py` — imbalance-aware stratified eval (AUPRC + bootstrap
  CI, precision/recall@K, HLA strata).
- `src/leakage.py`, `src/prime_leakage.py` — train/test overlap removal + paired
  bootstrap.

## Reproduce
Third-party tools/data are not redistributed. (1) Clone BigMHC, PRIME,
MixMHCpred, NeoRanking, ITSNdb. (2) Download BigMHC `im_train.csv` (Mendeley
dvmz6pkzvb) and the Müller neo-peptide matrix (figshare 147e67dde683fb769908).
(3) Run `python3 src/canonicalize.py` to build `data/clean/`, then the
eval/leakage scripts.

## Data sources & licenses
ITSNdb (Carri et al., Front. Immunol. 2023); BigMHC (Albert et al., Nat. Mach.
Intell. 2023, CC BY 4.0); PRIME/MixMHCpred (Gfeller lab); Müller/NeoRanking
(Müller et al., Immunity 2023, CC BY-NC); NEPdb (check source license).

## Related work (supersedes the framing here)
- Kim et al., T-SCAPE, *Science Advances* 2025 — leakage-controlled
  immunogenicity benchmark, explicit BigMHC leakage critique, allele-stratified
  eval.
- MHLAPre, *Brief. Bioinform.* 2024 — immunogenicity prediction fails off
  high-frequency training alleles.
- NeoTImmuML 2025; NeoMUST 2024 — leakage-controlled / independent-test eval.
