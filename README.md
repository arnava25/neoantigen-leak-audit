# NeoForge

Leakage controlled re-evaluation of HLA class I neoantigen immunogenicity
predictors (CD8+, single substitution neoepitopes).

**Status: exploratory / negative result, not a manuscript.** The central finding
(that leading immunogenicity tools are inflated by training data leakage on the
standard benchmark and generalize poorly off HLA-A\*02:01) was independently
published in stronger form while this was in progress; see *Related work*. This
repo is preserved for the reproducible pipeline and the data quality notes, not
as a novel claim.

## Key findings (ITSNdb, n=199)

**Benchmark contamination.** About 48% of ITSNdb test peptides (57% of positives)
appear in BigMHC's immunogenicity training set, and PRIME shows comparable
overlap against the shared public training pool. ITSNdb is clean as a *dataset*
but contaminated as a *benchmark* for tools trained on the same literature.

**Leakage free collapse off A\*02:01.** With leaked peptides removed, both tools
hold on A\*02:01 (AUROC about 0.78 to 0.84) but fall to chance or below on
non A\*02:01 (AUROC about 0.35 to 0.44). The apparent cross allele competence was
largely memorization.

**Cross allele miscalibration is tool specific.** Per allele recalibration
against a presented peptide background significantly improves BigMHC's
non A\*02:01 ranking (paired AUPRC delta +0.167) but degrades PRIME's. There is
no universal post hoc fix, and raw immunogenicity scores should not be pooled
across alleles for ranking.

## Caveats

Leakage free non A\*02:01 has only about 35 positives (underpowered, wide CIs).
PRIME's leakage set is approximated by BigMHC's training pool, which is
optimistic for PRIME. Single benchmark, two tools.

## Reusable

`src/canonicalize.py` does peptide and HLA canonicalization (mhcgnomes) plus the
data quality filters (the minimal epitope "unit trap" and the positive leakage
dedup).

`src/eval_baseline.py` does imbalance aware stratified evaluation (AUPRC with
bootstrap CI, precision and recall at K, HLA strata).

`src/leakage.py` and `src/prime_leakage.py` do train/test overlap removal plus a
paired bootstrap.

## Reproduce

Third party tools and data are not redistributed here. First, clone BigMHC,
PRIME, MixMHCpred, NeoRanking, and ITSNdb. Second, download BigMHC `im_train.csv`
(Mendeley dvmz6pkzvb) and the Muller neo peptide matrix
(figshare 147e67dde683fb769908). Third, run `python3 src/canonicalize.py` to
build `data/clean/`, then run the eval and leakage scripts.

## Data sources and licenses

ITSNdb (Carri et al., Front. Immunol. 2023). BigMHC (Albert et al., Nat. Mach.
Intell. 2023, CC BY 4.0, Mendeley dvmz6pkzvb). PRIME and MixMHCpred (Gfeller
lab). Muller / NeoRanking (Muller et al., Immunity 2023, CC BY-NC). NEPdb (check
source license before redistribution).

## Related work (supersedes the framing here)

Kim et al., T-SCAPE, *Science Advances* 2025: leakage controlled immunogenicity
benchmark, explicit BigMHC leakage critique, allele stratified evaluation.

MHLAPre, *Brief. Bioinform.* 2024: immunogenicity prediction fails off high
frequency training alleles.

NeoTImmuML 2025 and NeoMUST 2024: leakage controlled and independent test
evaluation.
