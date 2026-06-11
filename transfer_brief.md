# NeoForge Transfer Brief

Carrying machinery from evolveAGENT (de novo AMP generation) to NeoForge (discriminative ranking of HLA-I-presented neoantigens by CD8+ immunogenicity).

Each point is tagged:
- **[Grounded]** = specifically discussed/observed in this project's history
- **[Inferred]** = my generalization, not directly stated
- **[Not covered]** = history does not address it; not filling it in

A structural note that governs the rest: evolveAGENT was a generative search optimizing a population. NeoForge is a discriminative model scoring a fixed input list. The predictor, calibration, negative-set, and reference-set lessons transfer almost directly. The evolutionary / quality-diversity machinery transfers only if NeoForge later adds candidate generation or design, not for pure ranking. This distinction is flagged below where it changes the lesson.

---

## 1. Novelty / curiosity machinery

**How novelty was operationalized.** [Grounded] Novelty was a population-level proxy: `novelty = 1.0 - avg_sim`, where `avg_sim` was the average pairwise sequence similarity. The similarity metric was **Jaccard similarity over k-mer sets** of each sequence (precomputed k-mer sets, set intersection over union). The reference it was scored against was the **run's own archive** of recently seen sequences (`archive_kmers_ref`), with an explicit later proposal to instead score against the full APD database.

**Representation space.** [Grounded] Discrete k-mer sets of the raw sequence. There was no learned-embedding distance (no PWM-latent or autoencoder-space novelty). [Not covered] Any embedding-space novelty metric, because it was never used.

**Normalization.** [Grounded] Novelty was clamped to [0,1] (`max(0, min(1, avg_novelty))`). A separate diversity signal was population Shannon entropy `mean_H` normalized by `log2(20)` (max entropy over 20 amino acids). These were distinct inputs to the controller, not the same number.

**Tradeoff against other objectives.** [Grounded] Novelty was one of six motive weights (`amp`, `safety`, `stability`, `realism`, `novelty`, `curiosity`) that summed to 1.0 and were renormalized every generation by a cognitive controller (`AgentController`). It was bounded by a per-motive floor (`MIN_W`), and the related `curiosity` motive had a hard cap (`MAX_CUR`) to stop runaway exploration. There was also an `amp + safety >= 0.30` floor so the functional objectives could never be starved by exploration pressure. Surprise was defined as the prediction error between the ML AMP estimate and a heuristic AMP estimate; curiosity as novelty-weighted uncertainty; attentional salience as a memory that biased mutation toward high-fitness motifs. [Grounded] Critically, these mostly acted on **which offspring got generated** (mutation/operator selection), not on selection of survivors.

**Pathologies and fixes.**

- **Deceptive novelty from the wrong reference set.** [Grounded] This is the single most important novelty lesson. With novelty scored only against the run's own archive, the search kept converging on known-AMP-like scaffolds because that is what scored well on the AMP objective; novelty against the archive only pushed away from recent history, not away from the space you actually wanted distance from. The diagnosis at the time was explicit: "you are not discovering genuinely new chemistry, you are rediscovering known AMP patterns." Wang independently flagged ~40% APD similarity in outputs. The proposed fix was to compute Jaccard novelty against the **entire APD database** (a fixed external reference), not the moving archive, so the pressure pushed into genuinely uncharted space. [Inferred] In Lehman/Stanley terms this is classic deceptive/relative novelty, but the framing in-project was practical, not theoretical.

- **Novelty silently doing "less work."** [Grounded] Under NSGA-II selection, novelty was not one of the ranked objectives, so it only influenced mutation, not who survived. Diversity at selection time came from crowding distance and from MAP-Elites, not from the novelty term. The takeaway stated at the time: if you want novelty to actually shape outputs, it has to enter the selection/ranking step, not just operator choice.

- **Reward-hacking / sawtooth instability.** [Grounded] After an abandonment / niche-clearing event, fitness spiked without genuine improvement. Root cause was a penalty arithmetic artifact: crowding and historical-peak penalties reset on abandonment, so old sequences looked artificially better under the new penalty regime. Fixes: clear the MAP-Elites grid at abandonment, physically replace the population with 5% elites + 95% fresh randoms, and a `_skip_repro` flag skipping reproduction on the abandonment generation. K-mer Jaccard between the two resulting basins (mean ~0.076) was used afterward to confirm they were genuinely distinct, not the same peak relabeled.

- **Convergence collapse in long runs.** [Grounded] Single-scaffold-family convergence was a recurring risk; countermeasures were the multi-island architecture, niche clearing, and the MAP-Elites grid.

**Transfer note.** [Inferred] In NeoForge the role of novelty changes fundamentally. In evolveAGENT it was an *exploration reward* steering a search. Your bet is that dissimilarity-to-self is the immunogenicity signal, which makes it a *feature in a discriminative model*, not a search reward. The grounded lesson that survives the change is the reference-set lesson: distance only means what its reference set means. "Foreignness" must be computed against a fixed, comprehensive self-proteome reference (e.g. the human proteome as a k-mer/peptide set), never against the candidate list or any moving pool, or you will reward trivial intra-batch dissimilarity that has nothing to do with self-foreignness.

---

## 2. Multi-objective dynamics

**Objectives.** [Grounded] AMP activity, toxicity/safety, proteolytic stability, structural realism/naturalness, predicted E. coli MIC (a log10 regression model, used as a low-weight bonus, ~0.05), and novelty. Manuscript framing was "5 objectives (+ MIC + APD novelty)."

**What conflicted.** [Grounded] Stability was the persistent loser: it plateaued around ~0.65 while AMP score and safety both converged to ~0.83 to 0.84 and tracked each other closely. Stability behaved like a hard binding constraint the operator mix could not overcome, attributed to a small/weak stability training set. [Grounded] AMP-vs-novelty was the other inherent tension (maximizing predicted activity pulls toward known scaffolds; novelty pulls away), which is what produced the rediscovery problem in section 1.

**How they were combined.** [Grounded] Three approaches were used across the project's life: (a) a weighted composite fitness driven by the adaptive motive weights; (b) NSGA-II with crowding distance (Pareto); and (c) a MAP-Elites quality-diversity grid over charge x hydrophobicity, with Pareto-within-cells (up to 5 non-dominated sequences per cell). MAP-Elites with Pareto-within-cells was settled on as the right combination.

**Tradeoffs / what surprised.**
- [Grounded] MAP-Elites outperformed Pareto-only on max fitness (~0.87 vs ~0.78).
- [Grounded] The biggest surprise: the compound multi-objective fitness **masked a broken AMP predictor**. Early runs felt like they had a good gradient, but that gradient was actually being produced by stability, naturalness, toxicity, and realism pulling in different directions. The AMP model itself was outputting near-binary scores the whole time. The multi-objective structure hid the defect. This is the cross-cutting warning for NeoForge: a composite metric can look healthy while a component is dead.
- [Grounded] Fitness saturation/ceiling effects were a recurring drag (see section 5).

[Not covered] The exact static weight values of the final composite are not reliably in history; weights were adaptive (controller-set each generation), so there is no single fixed vector to port.

---

## 3. Predictor lessons

**Featurization that mattered.** [Grounded] Sequence input padded to a fixed `max_len` (20 for the 25-mer model) into a CNN. The decisive featurization choice was not architecture, it was the label and negative-set design (sections covered below). [Not covered] Specific layer/filter architecture details and any explicit calibration method (temperature scaling, etc.).

**Calibration behavior (the central lesson).** [Grounded] A retrained 13-mer AMP model collapsed to near-binary output (~0.02 or ~0.98). The exact phrasing at the time: "the fitness landscape has a cliff, not a gradient." A well-calibrated model outputs 0.3 / 0.5 / 0.7 and gives the search something to climb; the binary model gave a jump in gen 1-2 and then everything sat at the ceiling. The 25-mer model had good gradient by contrast. This is a calibration story, and for a ranking task it is existential: a near-binary scorer cannot rank.

**Where the CNN failed / was overconfident.** [Grounded] A model trained with shuffled-AMP negatives scored a poly-A sequence (`AAAAAAAAAAAAA`) at 0.927. It had learned a boundary that did not reject superficially AMP-like junk (poly-A, poly-L, poly-K). That is the predictor being reward-hackable: a sequence trivially distinct from real peptides scored as a strong hit.

**Root causes of the bad 13-mer model.** [Grounded] ~8,358 unique sequences, all from one run and one scaffold family, binary labels at a 0.5 threshold, and no hard negatives. The 25-mer model worked because it saw ~51k sequences from real databases (APD/DRAMP/DBAASP) with genuine positives and negatives across lengths.

**Evaluation mistake that fooled, and how it was caught.** [Grounded] Two of them. First, the composite-fitness masking described in section 2: the fix was to inspect the AMP model's raw outputs directly and to probe it with junk like poly-A, rather than trusting the composite. Second, the **circularity** problem: the internal CNN was trained on APD/DRAMP/DBAASP, and the external validators used to "confirm" results (CAMPR3, PeptideRanker, Veltri) draw on overlapping data and features, so high external agreement was partly tautological. This was raised by reviewers and by expert collaborators (van Hoek, de la Fuente), and was acknowledged in the manuscript. The resolution was honest framing plus treating wet-lab MIC/hemolysis as the only non-circular ground truth.

**Train/test leakage.** [Not covered] in the strict held-out-split-bug sense. [Inferred] The functional analogues that did occur: training the 13-mer model on one run's own outputs and then operating in that same narrow distribution (distributional overfit), and the circularity above. Treat these as the leakage lessons; there was no documented "I had a split bug that inflated accuracy" incident.

---

## 4. Data pipeline

**The negative set was the decision that mattered most.** [Grounded] Shuffled-AMP negatives poisoned the classifier (the poly-A 0.927 failure). The fix was real negatives: non-AMP bioactive peptides from UniProt (hormones, signal peptides, neuropeptides, real peptides that are definitively not antimicrobial) and/or experimentally inactive peptides (DBAASP MIC >= 128). The principle stated at the time: "negatives are where everyone struggles," and both classes must be real peptides so the model learns features that actually correlate with activity, not "is this a peptide at all."

**Continuous beats binary labels.** [Grounded] The AMP model was rebuilt (v4) on continuous MIC potency data rather than binary AMP/non-AMP, specifically to restore gradient and calibration.

**Do not mix a clean dataset with a poisoned one.** [Grounded] Combining the cleaned dataset with the old shuffled-negative dataset (`amp_dataset_v3`) would reintroduce the bad negatives; the cleaner set was kept on its own.

**Label-threshold choice is a quantity/purity tradeoff.** [Grounded] A strict MIC threshold (active <= 8, inactive >= 128) gave clean labels but only 1,722 sequences, too few for a CNN and prone to overfit. The threshold had to be relaxed to get enough data. There is a real tension between label purity and dataset size.

**Provenance and dedup hygiene.** [Grounded] Backup CSVs contained ~121k overlapping duplicate rows and had schema differences across run versions; a naive merge produced a corrupt combined file. The decision was not to merge but to archive separately. Counting peptides correctly across files mattered for the manuscript's headline number.

**Initialization as a data decision.** [Grounded] APD-seeded initialization produced ~40% APD similarity in outputs (a novelty/IP problem flagged by Wang); random initialization was preferred for novelty.

---

## 5. Dead ends (things deleted or abandoned, and why)

- **The 13-mer-specific model retrain.** [Grounded] Switching from 25-mer models to a 13-mer-specific retrain broke the system (binary collapse, fitness had to be redone). The eventual decision was to revert to the 25-mer models scoring 13-mers (accepting a length mismatch in exchange for a usable gradient), and ultimately to return to 25-mer work entirely. A model calibrated to a narrow length on its own run's data was worse than a broader model used slightly out of domain.

- **Shuffled-AMP negatives (`amp_dataset_v3`).** [Grounded] Replaced; see section 4.

- **Logistic compression as the fitness squashing function.** [Grounded] Created a hard ceiling that flattened the top of the scoring surface, preventing discrimination between good and exceptional sequences, and made early runs report fake ~1.0 best fitness. Loosened (beta/bias retuned), bonus capped at 1.15, softplus center raised, junk_cut raised, hard AMP gate (fitness 0 below 0.85 AMP score) added. The honest-numbers lesson: a ceiling makes everything look solved.

- **Pareto-only (NSGA-II without QD).** [Grounded] Underperformed MAP-Elites and let novelty fall out of selection entirely. Superseded by MAP-Elites + Pareto-within-cells.

- **The strict-MIC-threshold dataset (1,722 rows).** [Grounded] Too small; abandoned in favor of a relaxed threshold.

- **The naive CSV super-merge (`master_evolution_history_full_database.csv`).** [Grounded] Corrupt due to schema mismatch; flagged for deletion.

- **The experimental-data AMP model via `train_from_experimental.py`.** [Grounded] Produced the poly-A 0.927 result and was set aside; notably the failure was not fully diagnosed at the time, which is itself a lesson about moving on from a red flag without root-causing it.

[Not covered] Any deleted experiments beyond those above.

---

## Highest-value lessons most likely to transfer to immunogenicity

1. **The reference set is the entire bet.** [Grounded core, inferred application] Dissimilarity only means what you measure it against. In AMP, novelty-against-own-archive produced deceptive novelty and endless rediscovery of known scaffolds; the fix was a fixed external reference (APD). For NeoForge, "foreignness" must be distance to a fixed, comprehensive self-proteome (human proteome as a k-mer or peptide-window set), computed identically for every candidate, never against the candidate batch. If you get this right the central hypothesis is at least well-posed; if you get it wrong you will reward batch-relative novelty that is unrelated to self-foreignness. Corollary: do not let the same self-reference define both the foreignness feature and the immunogenicity label, or you rebuild the circularity problem in a new domain.

2. **Calibration is non-negotiable for a ranker.** [Grounded] A near-binary predictor gives a cliff, not a gradient, and a cliff cannot rank. Immunogenicity prediction is intrinsically a ranking problem, so you need genuine graded probabilities. Prefer continuous/graded labels (T-cell response magnitude, dilution-series readouts) over binary responder/non-responder where the data allows, and check calibration on the raw predictor directly, not through any composite score.

3. **Negatives, especially hard negatives, are the whole game.** [Grounded core, inferred application] Easy negatives gave a model that scored poly-A at 0.927; it learned "is this a peptide" instead of "is this active." The immunogenicity analogue is harsher: your hard negatives are peptides that *are* HLA-I-presented but *do not* elicit a CD8+ response, plus self peptides. They look almost identical to positives. If you train against random or non-presented peptides, you will learn presentation, which is already given to you as input, not immunogenicity. Source experimentally confirmed presented-but-non-immunogenic peptides as negatives.

4. **A composite or downstream metric can hide a dead component, and shared-data agreement is circular.** [Grounded] The multi-objective fitness masked a broken AMP model for a long time, and classifier agreement was partly tautological because internal and external models shared training data. For NeoForge, evaluate the immunogenicity predictor in isolation on held-out, experimentally labeled data, and keep the foreignness/self-distance feature computationally independent of the predictor's training features so you can measure whether foreignness adds real signal rather than re-deriving something the model already encodes.

5. **Quality-diversity beats a single ranked list if you care about coverage, and there is a real purity-vs-quantity tradeoff in the data.** [Grounded core, inferred application] MAP-Elites beat Pareto-only and prevented scaffold collapse; this transfers only if NeoForge surfaces a diverse set (across HLA alleles, anchor-residue patterns, epitope families) rather than scoring one fixed list, in which case a QD grid over meaningful descriptors beats a flat top-N that clusters on near-duplicates. And expect the same label-purity vs dataset-size tension that shrank the AMP dataset to 1,722 rows; immunogenicity datasets are small and noisy, so the threshold you set on "responder" will trade cleanliness against having enough data to train at all.