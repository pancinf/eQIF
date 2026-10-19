# eQIF: eQTL prediction with isolation forests

eQIF is an unsupervised machine-learning approach that combines 20 Enformer-derived principal-component scores using isolation forests to prioritize common variants by their predicted gene-regulatory impact. This repository contains the Snakemake preprocessing workflow, the R analysis notebook, and scripts for validation and regulatory analyses.

## Software requirements

**Environment to reproduce our results:** R **4.5.3** (2026-03-11), Snakemake **9.22.0**, bcftools **1.17**, CUDA **12.0** (build `cuda_12.0.r12.0`).

**R packages** (versions used):

| CRAN / other R packages | Version | Bioconductor packages | Version |
|---|---|---|---|
| `data.table` | 1.17.8 | `GenomicRanges` | 1.62.1 |
| `reshape2` | 1.4.5 | `rtracklayer` | 1.70.1 |
| `isotree` | 0.6.1.5 | `Rsamtools` | 2.26.0 |
| `dplyr` | 1.2.1 | `universalmotif` | 1.28.0 |
| `tidyr` | 1.3.2 | `TFBSTools` | 1.48.0 |
| `precrec` | 0.14.5 | `GenomicFeatures` | 1.62.0 |
| `ggplot2` | 4.0.3 | | |
| `ggpattern` | 1.3.1 | | |
| `gridExtra` | 2.3 | | |
| `ggseqlogo` | 0.2.2 | | |
| `patchwork` | 1.3.2 | | |
| `stringr` | 1.6.0 | | |
| `ggrepel` | 0.9.8 | | |
| `scales` | 1.4.0 | | |
| `RcppCNPy` | 0.2.15 | | |
| `purrr` | 1.2.2 | | |
| `optparse` | 1.8.2 | | |

Install R packages from CRAN (`install.packages()`) and Bioconductor (`BiocManager::install()`). Running `eQIF.ipynb` also requires **Jupyter and an R kernel (`IRkernel`)**.

**Additional command-line dependencies:** `samtools`, `tabix`, `bedtools` (v2.31.1), UCSC `liftOver` and `bigWigToBedGraph` ,and MEME Suite (`ame`, `fimo`) (v. 5.5.5). Install the **Sei framework** (https://github.com/FunctionLab/sei-framework) and its dependencies separately to generate Sei predictions; CUDA/GPU support is needed for the included Sei script.

## Repository layout

All commands and relative paths assume execution from `scripts/`:

```text
eQIF/
├── README.md
├── scripts/
│   ├── Snakefile
│   ├── eQIF.ipynb
│   ├── filterGTEx.R, filterGwas.R, getCds.R, liftEnformer.R
│   ├── runSeiAll.sh, ameEnrichment.sh, fimo.sh
│   ├── liftOver, bigWigToBedGraph   # UCSC executables
│   └── sei-framework/              # external Sei installation
└── data/
    ├── raw/                        # source data (below)
    └── outputs/                    # generated data and results
```

## Required input data

Place the following files **under `data/raw/`** before execution (paths and filenames are fixed in the scripts unless otherwise noted):

| Analysis | Required paths within `data/raw/` |
|---|---|
| Variant preprocessing | `gnomAD/<input>.vcf.gz` (path set via Snakemake config), `reg/hla.bed`, `reg/gencode.v44.annotation.gtf` |
| Fine-mapped eQTL / GWAS | `gtexEqtl/CAVIAR_Results_v8_GTEx_LD_ALL_NOCUTOFF_with_Allele.txt.gz`, `gwasFinemap/v2.1/credible_set.txt` |
| Enformer / coordinate conversion | `enf/*.tsv.bgz` (Enformer scores, including PC0–PC19), `chain/hg19ToHg38.over.chain`, `chain/hg19ToHg38.over.chain.gz` |
| Other prediction scores | `linsight/LINSIGHT.bw`, `dvar/hg19_DVAR.score.gz`, `cadd/whole_genome_SNVs_38.tsv.gz`, `avi/alphagenome_variant_impact_score_snvs.tsv.gz` |
| Conservation / regulatory regions | `phylos/hg38.phyloP100way.bw`, `phylos/hg38.phyloP447way.bw`, `phylos/hg38.phyloP447wayPrimates.bw`, `reg/GRCh38-cCREs.bed`, `reg/consensus_footprints_and_collapsed_motifs_hg38.bed` |
| MPRA / locus examples | `mpra/GTExMpra.txt`, `candidateLoci/mx2.txt`, `candidateLoci/sort1.txt` |
| Motif analyses | `ref/hg38.fa`, `motifs/MA0139.2.meme`, `motifs/MA0098.4.meme`, `ame/clusters.tab`, `ame/annotation_table.txt` |

Note that the file links can be found in our publication.

**Additional generated inputs for the notebook:** 
Sei sequence-class scores for chromosomes 1–22 at `data/outputs/Sei/gn_v3_red_inf_id2_seiAll_chr{1..22}.sequence_class_scores.npy`, plus CADD/AVI lookup files created by the notebook or provided in advance (`data/outputs/cadd/cadd_hits.tsv`, `data/outputs/avi/avi_hits.tsv`). Sei score rows must match the gnomAD variant order.

## Running the workflow

1. Download the source data to the paths above, install the required software, and create output directories as needed (including `data/outputs/eval/ame/` and `data/outputs/eval/fimo/`).
2. From the `scripts/` directory, run preprocessing:

   ```bash
   snakemake -s Snakefile --cores 3 \
     --config gnomadFilePathArg="../data/raw/gnomAD/<input>.vcf.gz" \
              bcftoolsPathArg="bcftools"
   ```

3. Generate Sei scores if not already available (requires `scripts/sei-framework/`):

   ```bash
   bash runSeiAll.sh
   ```

   The script sets `CUDA_VISIBLE_DEVICES=1`; adjust it for your GPU configuration.

4. Open `eQIF.ipynb` with the Jupyter R kernel **from `scripts/`** and execute the cells in order. Set `runFirst <- TRUE` in the CADD/AVI annotation cell on the first run to generate the tabix lookup files (the default is `FALSE`). The notebook runs the isolation-forest analyses and downstream evaluations, including AME/FIMO.

**Outputs:** Preprocessing files are written to `data/outputs/`; figures, summary tables and exported eQIF scores are written to `data/outputs/eval/`.
