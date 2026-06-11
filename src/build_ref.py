import pandas as pd
test = pd.read_csv("data/clean/itsndb_test.csv")
mul  = pd.read_csv("data/clean/muller_train.csv")
neg  = mul[mul.label == 0]

N = 200
test_alleles = test.hla.unique()
ref = []
for a in test_alleles:
    pool = neg[neg.hla == a]
    if len(pool):
        ref.append(pool.sample(min(N, len(pool)), random_state=0))
ref = pd.concat(ref, ignore_index=True)
ref[["hla", "peptide"]].rename(columns={"hla": "mhc", "peptide": "pep"}).assign(tgt=0)\
   .to_csv("data/clean/bigmhc_ref_input.csv", index=False)

counts = neg.groupby("hla").size()
have = [a for a in test_alleles if counts.get(a, 0) >= 20]
nonA2_have = [a for a in have if a != "HLA-A*02:01"]
print(f"reference peptides: {len(ref)} across {ref.hla.nunique()} alleles")
print(f"test alleles with >=20 Müller refs: {len(have)}/{len(test_alleles)} "
      f"(non-A*02:01: {len(nonA2_have)})")
print("non-A*02:01 test alleles WITHOUT enough refs:",
      [a for a in test_alleles if a != "HLA-A*02:01" and counts.get(a, 0) < 20])
