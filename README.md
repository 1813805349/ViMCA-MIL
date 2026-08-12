# ViMCA-MIL

**Vi**ral **M**utation–**C**linical **A**ssociation analysis based on a gated-attention **M**ultiple **I**nstance **L**earning model
<img width="2126" height="616" alt="image" src="https://github.com/user-attachments/assets/09cec173-c381-4bc5-80ae-653dba06e56b" />
</br>

## What is ViMCA-MIL?

ViMCA-MIL is an integrated analysis framework that screens for SARS-CoV-2 mutations
which simultaneously confer viral evolutionary advantages and impact host clinical
phenotypes. It integrates large-scale genomic surveillance data with matched
virus–patient clinical records, and prioritizes candidate mutations through a
gated-attention multiple-instance learning (MIL) model.

Briefly, the framework consists of four steps:

1. **Evolutionary analysis.** Variant and mutation frequency trajectories are
   characterized across major SARS-CoV-2 variants using ~15.5 million genomic
   sequences from GISAID (Dec 2019 – Jan 2026) together with locally sequenced
   genomes from West China Hospital, Sichuan University.
   Recurrent mutations are clustered into evolutionary clusters to distinguish
    *stable* from *sporadic* mutations.
3. **Fitness and transmission analysis.** Relative viral fitness (R/R<sub>A</sub>)
   is estimated with the [PyR0](https://github.com/broadinstitute/pyro) workflow on
   the UShER phylogeny, and transmission selection coefficients are inferred from
   regional genomic surveillance data following
   [Lee *et al.*](https://www.nature.com/articles/s41467-024-55593-0).
   EVEscape-derived fitness effects are additionally computed per protein.
4. **Clinical phenotype analysis.** For a cohort of 490 patients with matched viral
   genomes and clinical records (78 blood routine, blood chemistry, immune cell and
   cytokine/chemokine features), features are normalized (inverse normal
   transformation), adjusted for age, sex and a Charlson comorbidity index (CCI)
   weighted clinical burden score, and filtered by temporal correlation with the
   course of the pandemic, yielding 38 temporally correlated traits.
5. **Gated-attention MIL modeling.** Each patient is treated as a *bag* of the
   mutations detected in the matched viral genome; each mutation is an *instance*
   represented by 15 mutation-level features. The model predicts each patient-level
   clinical feature value from the attention-weighted aggregation of its mutations,
   and the attention scores are used to prioritize mutations without requiring
   mutation-level labels. Robustness is assessed with 100 independent resampling
   runs followed by 10 full-data refittings to obtain mean attention scores.

Using this framework, we identified high-frequency stable mutations in nonstructural
proteins (Nsps), notably **Nsp1:S135R** and **Nsp13:R392C** emerged in Omicron,
which are associated with serum TNF-α levels and were experimentally validated to
enhance viral replication while attenuating pulmonary pathogenicity.

</br>

## Contents

[Repository structure](#repository-structure)  
[Requirements](#requirements)  
[Input](#input)  
[Output](#output)  
[Tutorial](#tutorial)  
[Mutation features](#mutation-features)  
[Contact](#contact)  
[Citation](#citation)

</br>

## Repository structure

```
ViMCA-MIL/
├── README.md
└── script/
    ├── 01Evolutionary_analysis_variant_mutations_v2.Rmd   # Variant/mutation evolutionary dynamics
    ├── 02Fitness_analysis_v2.Rmd                          # Relative fitness & transmission selection
    ├── 03Clinical_phenotype_analysis_v2.Rmd               # Clinical feature preprocessing & temporal analysis
    ├── MIL_bootstrap_multi_trait_260714.py                # Gated-attention MIL model (single trait per run)
    ├── run_MIL.sh                                         # Batch launcher for all traits + summary merging
    └── count_nm_sm_by_sequence_gene.py                    # Per-sequence/per-gene Nm & Sm counting (dN/dS)
```

</br>

## Requirements

### Python (MIL model)

Developed with `Python 3.8+` and `PyTorch`. Required packages:

    torch (>= 1.10),
    numpy (>= 1.21),
    pandas (>= 1.3),
    scikit-learn (>= 1.0),
    scipy (>= 1.7),
    tqdm (>= 4.60)

```bash
pip install torch numpy pandas scikit-learn scipy tqdm
```

A CUDA-capable GPU is recommended; the code automatically falls back to CPU.

### R (evolutionary / fitness / clinical analysis)

Developed with `R 4.2.1`. Main packages used:

    data.table, openxlsx, readxl, dplyr, tidyr, purrr, stringr, lubridate,
    ggplot2, ggpubr, ggrepel, patchwork, cowplot, ggsci, RColorBrewer,
    ComplexHeatmap, circlize, eulerr, binom, caret, bestNormalize, broom,
    future.apply,  Biostrings, clusterProfiler

External tools invoked by the analysis notebooks:

- [fastp](https://github.com/OpenGene/fastp) (v0.23.2) — read preprocessing
- [bowtie2](https://github.com/BenLangmead/bowtie2) (v2.4.4) — read alignment
- [Nextclade](https://clades.nextstrain.org/) — mutation calling and lineage assignment
- [PyR0](https://github.com/broadinstitute/pyro) — relative fitness estimation
- [epi-covar](https://www.nature.com/articles/s41467-024-55593-0) — transmission selection coefficients
- [EVEscape](https://github.com/broadinstitute/EVEscape) + [EVcouplings](https://github.com/debbiemarkslab/EVcouplings) — fitness-effect scores

</br>

## Input

> **Note:** The scripts contain hard-coded working paths (e.g.
> `SCRIPT`, `TRAIT_FILE`, `OUT_DIR` in `run_MIL.sh`; `GENO_PATH`, `PHENO_PATH` in
> the MIL Python script). Please adapt these paths to your own environment before
> running.

#### 1. Mutation feature matrix — `geno_feature.csv`

A long-format data frame in which each row is one mutation detected in one patient,
containing at least the following columns:

| Column | Description |
| --- | --- |
| `sample_id` | Patient/sample identifier (defines the MIL *bag*) |
| `mutation_id` | Mutation identifier (e.g. `Spike:D614G`) |
| 15 feature columns | See [Mutation features](#mutation-features) |

#### 2. Clinical phenotype matrix — `nor_pheno_by_Age_Gender_CCI.csv`

A data frame with one row per patient (`sample_id`) and one column per clinical
feature. Values are normalized, covariate-adjusted (age, sex, CCI) and
standardized, as described in `03Clinical_phenotype_analysis_v2.Rmd`.

#### 3. Trait list — `traits.txt`

A plain-text file listing one clinical trait (column name of the phenotype matrix)
per line. `run_MIL.sh` launches one MIL job per trait.

#### 4. Multiple sequence alignment (optional, dN/dS analysis)

A FASTA (optionally gzipped) multiple sequence alignment of SARS-CoV-2 genomes
(e.g. the GISAID MSA), plus a gene/feature annotation table, used by
`count_nm_sm_by_sequence_gene.py`.

</br>

## Output

Running `run_MIL.sh` produces, for each trait, under `OUT_DIR`:

- `predictions/` — observed vs. predicted values for each resampling run
- `attention_each_run/` — mutation attention scores per run
- `model_weights/` — trained model weights
- `summary_each_trait/Summary_*.csv` — per-trait performance summary
  (mean Pearson *r* over the 100 resampling runs, etc.)
- `Summary_All_Traits.csv` — merged summary across all traits
- `logs/` — per-trait log files

The final mutation–clinical feature association is quantified by the **mean
attention score** of each mutation for each clinical feature obtained from the 10
full-data refittings. Mutations with attention scores above the mean across all
evaluated mutations are defined as high-attention mutations for that feature.

</br>

## Tutorial

### Step 1 — Evolutionary analysis of variants and mutations

`script/01Evolutionary_analysis_variant_mutations_v2.Rmd`

- Aggregates Pango lineages into major variant groups (e.g. `BA.2*`, `JN*`,
  `XBB.1.5*`, `XDV*`, …) and computes weekly variant frequencies (global vs. China).
- Retains 1,161 recurrent nonsynonymous mutations (frequency > 5% at any time point).
- Builds a mutation–variant frequency matrix and clusters mutations into 12
  evolutionary clusters using hierarchical clustering (`hclust`, `ward.D2`),
  distinguishing stable from sporadic mutations.
- Computes length-normalized mutation density per viral protein and plots mutation
  frequency trajectories.

### Step 2 — Viral fitness and transmission analysis

`script/02Fitness_analysis_v2.Rmd`

- **Relative fitness**: runs the PyR0 workflow (`preprocess_usher.py`,
  `mutrans.py`) on the UShER tree from the UCSC Genome Browser.
- **Transmission selection coefficient**: computes region-specific temporal
  prevalence trajectories (`epi-covar.py`) and infers selection coefficients
  (`epi-inf-parallel.py`).
- Maps ORF1a/ORF1b mutations to Nsp coordinates and integrates EVEscape-derived
  fitness effects.

### Step 3 — Selection pressure (dN/dS) analysis

`script/count_nm_sm_by_sequence_gene.py`

Counts nonsynonymous (Nm) and synonymous (Sm) substitutions per sequence and per
gene directly from a SARS-CoV-2 MSA, relative to the Wuhan-Hu-1 reference, which
are used to estimate dN/dS with the Nei–Gojobori method:

```bash
python count_nm_sm_by_sequence_gene.py \
    --msa gisaid_msa.fasta.gz \
    --features gene_annotation.csv \
    --out nm_sm_counts.csv \
    --skip-reference
```

### Step 4 — Clinical phenotype preprocessing

`script/03Clinical_phenotype_analysis_v2.Rmd`

- Computes the Charlson-weighted clinical burden score (modified CCI by Quan *et
  al.*) from inpatient diagnostic records.
- Applies inverse normal transformation (`bestNormalize`), stepwise covariate
  adjustment (age, sex, CCI; `caret`) and standardization to all 78 clinical
  features.
- Calculates Pearson temporal correlations between each clinical feature and the
  infection date; the 38 significant features (P < 0.05) are retained as the MIL
  phenotype set (16 severe-related, 14 non-severe-related, 8 unknown).

### Step 5 — Gated-attention MIL modeling

`script/MIL_bootstrap_multi_trait_260714.py` + `script/run_MIL.sh`

Single-trait usage:

```bash
CUDA_VISIBLE_DEVICES=0 python MIL_bootstrap_multi_trait_260714.py \
    --trait "TNF-a" \
    --output_dir ./multi_traits_results \
    --n_boot 100 \
    --n_full 10 \
    --epochs 150 \
    --base_seed 42 \
    --r_threshold 0 \
    --save_phase1_attention 1 \
    --save_phase1_pred 1 \
    --save_phase2_model 1 \
    --save_phase2_attention 1 \
    --save_phase2_pred 1
```

The pipeline runs in two phases:

- **Phase 1 — stability screening**: 100 independent resampling runs; in each run
  80% of patients are used for training and 20% for validation. Models are trained
  for 150 epochs with MSE loss and the AdamW optimizer (lr = 5×10⁻⁴, weight decay
  = 0.01).
- **Phase 2 — refitting**: for traits whose mean validation Pearson correlation
  exceeds the threshold (`--r_threshold`), the model is refit 10 times on all
  matched samples with different seeds; mean mutation attention scores are
  extracted (see [MIL model](#mil-model)).

To run all traits in parallel (6 jobs per batch by default), edit the paths at the
top of `run_MIL.sh` and run:

```bash
bash run_MIL.sh
```

Downstream, traits whose mean validation Pearson correlation exceeds the upper
bound of the 95% confidence interval across all evaluated traits are retained for
mutation prioritization (13 traits in our study, with serum TNF-α being the most
predictable).

</br>

## Mutation features

Each mutation instance is represented by 15 features:

| Category | Features | Source |
| --- | --- | --- |
| Conservation | `phastCons`, `phyloP` (max over the affected codon, 119 coronavirus genomes) | UCSC Genome Browser |
| Fitness & transmission | `relative_fitness` (R/R<sub>A</sub>), `selection_coefficient`, `evo_idx` (EVEscape fitness effect) | PyR0 / Lee *et al.* / EVEscape |
| Protein stability | `delta_DDG_Env`, `delta_DDG_Int` (predicted ΔΔG under environmental/intracellular conditions) | Cov2Var |
| Physicochemical properties | `Molecular_weight`, `Theoretical_PI`, `Extinction_coefficients`, `Aliphatic_index`, `grand_average_of_hydropathicity` | Cov2Var |
| Functional impact | `Protein_Func` (pathogenicity), `SIFT_Probability`, `PROVEAN_Score` | Cov2Var / SIFT / PROVEAN |

</br>

## Contact

For questions about the code, please open an issue or contact
Lu Chen (luchen@scu.edu.cn) or Kepan Linghu (lhkp5457@163.com).

</br>

## Citation

If you find ViMCA-MIL useful, please cite our paper:

> Wei H.-C.\*, Linghu K.\*, Yang H.\*, Yang Q.\*, Huang X.\*, Wang Y.-H., *et al.*
> Effects of high-frequency clinical mutations on SARS-CoV-2 replication and
> virulence. *(under submission)*

</br>

