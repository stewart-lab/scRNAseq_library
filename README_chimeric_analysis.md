# Human–pig chimeric scRNA-seq analysis

This pipeline analyzes 10x 3′ single-cell data from samples that mix **human and pig cells** (e.g. the ChimeraTLS 20% human / 80% pig somite cultures). It:

1. aligns every read to a **combined human + pig reference**, so each read is assigned to the species it fits best;
2. **calls each cell barcode human, pig, mixed (human+pig doublet) or nuclear doublet**;
3. measures **ambient (background) RNA** using the other species as a tracer, giving the SoupX contamination fraction ρ;
4. compares the species calls with STARsolo's own cell calls;
5. tests whether other-species reads in a cell are **more than ambient RNA and mismapping can explain**, which is the collaborator's question about human RNA being "taken up" by pig cells and vice versa.

| Step | Script | What it does |
|---|---|---|
| 0 | `extend_3prime_gtf.py` | Extends pig gene 3′ ends (the pig annotation lacks many 3′ UTRs) |
| 1 | `01_build_combined_reference.sh` (+ `prefix_gtf.py`) | Builds the combined human+pig FASTA, GTF and STAR index |
| 2 | `run_aligner.sh` + `config.json` | Runs STARsolo against the combined index: all lanes of the library in one run (`COMBINE_LANES: true`), or one run per lane |
| 2b | `02_merge_lanes.py` | Combines existing per-lane STARsolo outputs of one library without re-aligning |
| 3 | `03_classify_species.py` | Species call per barcode, ρ, nuclear-doublet call |
| 4 | `04_compare_filtered_barcodes.py` | Compares the species calls with STARsolo's filtered (cell-called) barcodes |
| 5 | `05_cross_species_excess.py` | Tests for other-species RNA beyond ambient RNA and mismapping |
| 6 | `run_aligner.sh` (single genome) → `run_seuratv5.sh` | Per species: whole-library alignment to one genome, SoupX with ρ from step 3, subset to that species' barcodes, then the standard Seurat pipeline |

---

## Step 0–1: Build the combined reference

**Location:** `/w5home/bmoore/genomes/human-pig_chimera/`

| File | Contents |
|---|---|
| `human/GRCh38.fa` | GRCh38 primary assembly **plus the reporter contig `GRCh38_MEOX1_reporter`** (NLuc-P2A-tdTomato + SV40 polyA, then an FRT-flanked NeoR cassette) |
| `human/GRCh38.108.gtf` | Ensembl 108 human annotation (filtered) + reporter genes `GRCh38_NLuc_P2A_tdTomato`, `GRCh38_NeoR` |
| `pig/Sscrofa11.fa` | **Sscrofa11.1 top-level** assembly: 613 sequences (chromosomes 1–18, X, Y, MT + scaffolds) |
| `pig/Sscrofa11_v1.109_3pext.gtf` | Ensembl 109 pig annotation with 3′ ends extended by `extend_3prime_gtf.py` |
| `combined.fa`, `combined.gtf` | Human + pig, concatenated |
| `star_solo/` | Combined STAR index: 71,329 genes (39,548 human + 31,781 pig) |

### Pig 3′ extension (step 0, optional)

Ensembl pig genes (and other non-model species) often stop at the stop codon with no 3′ UTR (about 18% of protein-coding genes, including FOXC2 and SOX2). 10x 3′ reads pile up past the annotated end and would be thrown away as intergenic. `extend_3prime_gtf.py` extends each gene's 3′ end using read coverage from 3′ scRNA-seq BAMs and/or NCBI RefSeq 3′ ends. It never extends into the next gene on the same strand. 

usage: `extend_3prime_gtf.py --gtf in.gtf --out out.gtf --report report.tsv \
                            [--refseq-gff refseq.gff] \
                            [--coverage-out cov.tsv.gz | --coverage-in cov.tsv.gz] \
                            [aligned.sam|bam ...]`

### What the build does

Run `bash 01_build_combined_reference.sh` (needs STAR on PATH, e.g. `conda activate star_env`, and about 60–100 GB RAM). Source files, all local:

| Input | File |
|---|---|
| Human FASTA | `human_genome_38/Homo_sapiens.GRCh38.dna_sm.primary_assembly_tdTomato_NanoLuc.fa` |
| Human GTF | `human_genome_38/Homo_sapiens.GRCh38.108.filtered.gtf` (includes the reporter genes) |
| Pig FASTA | `pig_genome/ensemble/Sus_scrofa.Sscrofa11.1.dna.primary_assembly.allchr.fa` |
| Pig GTF | `pig_genome/ensemble/Sus_scrofa.Sscrofa11.1.109.chr.filtered.3pext.gtf` (output of step 0) |

1. **Prefix every contig and gene with its species** (`GRCh38_`, `Sscrofa11_`) using `awk` for the FASTA and `prefix_gtf.py` for the GTF. Every contig, gene ID and gene name then shows its species, and nothing collides between the two genomes.
2. **No biotype filtering.** Human and pig annotations use different biotype vocabularies. Filtering to `protein_coding` would remove a different fraction of each species and bias the species ratio.
3. **Concatenate and index:**
   ```bash
   STAR --runMode genomeGenerate --runThreadN 24 --genomeDir star_solo \
        --genomeFastaFiles combined.fa --sjdbGTFfile combined.gtf \
        --sjdbOverhang 149 --genomeSAsparseD 2 --limitGenomeGenerateRAM 90000000000
   ```
   Needs about 40–60 GB RAM and a few hours (`--genomeSAsparseD 2` lowers the memory needed). Run it on a cluster node.



### Checks after building (both caught real problems)

- **Every pig exon must be indexed.** The FASTA and GTF must use the same chromosome names:
  ```bash
  grep -c 'not found in Genome fasta' star_solo/Log.out   # must be 0
  cut -f1 star_solo/geneInfo.tab | tail -n +2 | sed 's/_.*//' | sort | uniq -c
  ```
  An earlier build used the Berkshire contig-level pig FASTA (`Sscrofa11_LUXW…` contigs) with the chromosome-based Sscrofa11.1 GTF. STAR dropped all 552,412 pig exons, so pig reads aligned but were counted under no gene. The current build uses the Sscrofa11.1 top-level FASTA and indexes every gene.
- **No duplicate contig names** (the script checks this).

---

## Step 2: Align to the combined reference

**Script:** `source run_aligner.sh config.json` (interactive Docker wrapper around `pre_pipeline/STAR.sh`)

Settings in `config.json`:
```json
"GENOME_DIR": "/w5home/bmoore/genomes/human-pig_chimera/",
"GENOME_INDEX_DIR": "star_solo/",
"fastq_alignment": { "SAMPLE_1": { "NAME": "ChimeraTLS_20H80P_3pext", "NUM_LANES": 2, ... } }
```
At the prompts, answer *no* to Parse and *no* to building an index, and *yes* to "Align to GENOME_INDEX_DIR".

### Lanes: combine them (`COMBINE_LANES`)

The FASTQs `…_S3_L001_…` and `…_S3_L002_…` are **one 10x library sequenced on two flow-cell lanes, so both lanes contain the same cells**. Each lane reads a random half of the same barcoded cDNA. On the ChimeraTLS data:
- 18,063 pig barcodes are shared between lanes; about 93 would be expected by chance if the lanes held different cells;
- per-barcode UMI counts correlate at r = 0.998;
- **76,992 identical molecules** (cell barcode + UMI + gene) appear in both lanes, against 0 when barcodes are scrambled between lanes.

Treating the lanes as two samples counts every cell twice downstream (merged objects, Harmony across lanes, inflated DE statistics), and lets a per-lane `--min-umi` give the same cell different calls. Combine the lanes in one of two ways:

- **`"COMBINE_LANES": true`** in the sample's config entry (recommended for new alignments). `pre_pipeline/STAR.sh` and `STAR_genomedir.sh` then pass all lanes to **one** STARsolo run (comma-separated FASTQs) and write to `shared_mount/<NAME>/`. STARsolo collapses UMIs across lanes, so counts are exact. Without it, each lane is aligned separately into `shared_mount/<NAME>_lane<i>/` as before.
- **`02_merge_lanes.py`** for lanes that are already aligned. It sums the per-lane raw matrices, checking that genes and barcodes match, and takes the union of the lanes' STARsolo cell calls as `filtered/`. Output goes to `<NAME>/Solo.out/<feature>/{raw,filtered}`, the same layout as a `COMBINE_LANES` run. About 3 minutes per feature.
  ```bash
  python3 02_merge_lanes.py \
      --lanes shared_mount/ChimeraTLS_20H80P_mapped2pigonly_lane1 shared_mount/ChimeraTLS_20H80P_mapped2pigonly_lane2 \
      --features Gene GeneFull_Ex50pAS
  # -> shared_mount/ChimeraTLS_20H80P_mapped2pigonly/Solo.out/...
  ```
  Caveat: a molecule read in both lanes is counted once per lane. On this library that inflates counts by **3.2%, almost the same for every gene** (2.9% for the least-expressed genes, 3.2% for the most), so it doesn't affect species calls, ρ, SoupX, clustering or DE. UMI totals, and so `--min-umi`, read about 3% high.

For the Seurat step, set `"COMBINE_LANES": true` for the sample either way. `get_data.py` then reads `shared_mount/<NAME>/` and ignores per-lane folders; with `false` it reads only the per-lane folders. After combining, run 03, 04 and 05 **once** on the combined folder, with one `--min-umi`. Depth is roughly doubled, so pick the cutoff again with step 4.

STARsolo settings (10x 3′ v4):

| Setting | Value | Why |
|---|---|---|
| `--soloCBwhitelist` | `3M-3pgex-may-2023_TRU.txt` | 10x 3′ v4 barcode whitelist |
| `--soloFeatures` | `Gene GeneFull GeneFull_Ex50pAS SJ Velocyto` | Exon-only and intron-inclusive matrices, both used downstream |
| `--soloMultiMappers` | `EM` | Writes `UniqueAndMult-EM.mtx` in addition to unique-only `matrix.mtx` |
| `--soloCellFilter` | `EmptyDrops_CR` | STARsolo cell calling, used for comparison in step 4 |
| `--soloUMIdedup` / `--soloUMIfiltering` | `1MM_CR` / `MultiGeneUMI_CR` | Cell Ranger-like |
| `--outSAMattributes` | `NH HI AS nM CB UB GX GN` | Needed for BAM-level checks (barcodes, gene tags, mismatches) |

**Data location (moved after alignment):**
`/w5home/bmoore/scRNAseq/LiFangChu/10x_somite_data/chimera_data/ChimeraTLS_20H80P_3pext_lane{1,2}/`

**Which matrix to use where:**

| Use | Matrix |
|---|---|
| Species calls and cross-species tests | **`GeneFull_Ex50pAS/raw/matrix.mtx`**: intron-inclusive, **unique reads only** |
| Exon-only counts for the intronic check | `Gene/raw/matrix.mtx` (same run) |
| Gene counts within a species | `UniqueAndMult-EM.mtx`, or the single-species alignments (see next steps) |

Do **not** use the EM matrix for cross-species work. EM splits ambiguous human/pig reads in proportion to each cell's own species, which adds uncertain cross-species counts. Use the *raw* matrices, not *filtered*: the empty droplets are needed to measure ambient RNA.

---

## Step 3: Classify species (`03_classify_species.py`)

```bash
D=/w5home/bmoore/scRNAseq/LiFangChu/10x_somite_data/chimera_data/ChimeraTLS_20H80P_3pext_lane1
python3 03_classify_species.py \
    --solo-dir      $D/Solo.out/GeneFull_Ex50pAS/raw \
    --exon-solo-dir $D/Solo.out/Gene/raw \
    --gtf           /w5home/bmoore/genomes/human-pig_chimera/combined.gtf \
    --human-purity 0.8 --pig-purity 0.9 --min-reporter-umi 2 \
    --outdir        $D/classified_cells_h0.8_p0.9_nd
```
(`--gtf` is only needed for the mitochondrial check; without it `mt_call_concordance` is reported as null.)

### Species call
For every barcode, UMIs are summed into human and pig totals using the gene prefixes. Barcodes with at least `--min-umi` (500) species-assignable UMIs are called:

| Call | Rule |
|---|---|
| `human` | human fraction ≥ `--human-purity` |
| `pig` | pig fraction ≥ `--pig-purity` |
| `mixed` | neither: mostly human+pig doublets |
| `unassigned` | fewer than `--min-umi` UMIs: empty droplets |

**Why separate cutoffs per species:** ambient RNA comes mostly from the majority species (the soup here is about 93% pig). A clean human cell therefore sits at about **0.89 human**, while a clean pig cell sits above 0.98 pig. A symmetric 0.9 cutoff discarded about 80% of human cells into `mixed`. Use **human ≥ 0.8, pig ≥ 0.9**, and check the histogram in `barnyard.png`. The script warns when a cutoff is at or above where clean cells are expected.

### Reporter override
`GRCh38_NLuc_P2A_tdTomato` and `GRCh38_NeoR` exist only in the human cells. A barcode with ≥ `--min-reporter-umi` reporter UMIs is called human regardless of purity. In practice this changes almost nothing (0–1 cells per lane):
- tdTomato appears to be knocked into the **MEOX1** locus. It's detected in only about 6% of human cells, and its detection follows MEOX1 expression.
- Its reads pile up correctly just before the SV40 polyA, inside the annotated region, so the low detection is not an annotation problem.
- NeoR has essentially zero reads, probably because it was excised by Flp.

### Ambient RNA and ρ for SoupX
The soup is barcodes with 1–100 UMIs (`--soup-max-umi`, SoupX's default range). The other species acts as a tracer:

```
rho_pig   = (human UMIs / all UMIs, pooled over pig cells)   / human share of soup
rho_human = (pig UMIs   / all UMIs, pooled over human cells) / pig share of soup
```

ρ is the **total** ambient fraction, including the cell's own species. The other species just makes it visible. Pass `rho_pig` and `rho_human` to SoupX `setContaminationFraction()` on the single-species alignments. Ambient RNA grows sub-linearly with cell size (about N^0.6), so smaller cells are more contaminated and ρ is an average. A noisy per-cell estimate is in the `rho_cell` column.

### Nuclear doublets (needs `--exon-solo-dir`)
Some barcodes pass the purity cutoff but carry a **piece of the other species' nucleus**. A cell is called `nuclear_doublet` when:
1. it has more foreign UMIs than ambient RNA predicts for its size (negative-binomial model, log μ = a + b·log N, BH q < `--excess-fdr`); **and**
2. those foreign reads are **intronic**: the exon-only/intron-inclusive ratio is significantly below that of ordinary cells (binomial, BH q < `--intronic-fdr`).

Ambient and mismapped reads are spliced (about 89% exonic). Nuclear material brings unspliced pre-mRNA. These cells get `call = nuclear_doublet`, are written to `barcodes_nuclear_doublet.txt`, are **removed from `barcodes_pig.txt` / `barcodes_human.txt`**, and are left out of ρ. `call_purity` keeps their original call.

### Outputs
| File | Contents |
|---|---|
| `species_assignment.csv.gz` | One row per barcode: species UMI sums, `frac_human`, `call_purity`, `call`, `reporter_override`, `rho_cell`, `foreign_expected`, `foreign_excess_q`, `foreign_exonic_frac`, `intronic_q`, `nuclear_doublet` |
| `barcodes_{human,pig,mixed,nuclear_doublet}.txt` | Barcode lists for subsetting |
| `summary.json` | Counts, purity cutoffs, ambient fractions, ρ, doublet-rate estimate, nuclear-doublet statistics |
| `barnyard.png` | Human vs pig UMIs per barcode + histogram of human fraction, with the cutoffs drawn |

---

## Step 4: Compare with STARsolo cell calls (`04_compare_filtered_barcodes.py`)

```bash
python3 04_compare_filtered_barcodes.py \
    --classified-dir $D/classified_cells_h0.8_p0.9_nd \
    --filtered-dir   $D/Solo.out/GeneFull_Ex50pAS/filtered \
    --outdir         $D/classified_cells_h0.8_p0.9_nd/filtered_comparison
```

This checks how the species calls overlap with STARsolo's EmptyDrops cell calls. Outputs:
- `overlap_summary.tsv`: per call, the number and percentage found in the filtered set, with UMI medians
- `filtered_not_called.tsv`/`.txt`: filtered barcodes the species call left unassigned
- `called_not_filtered_<call>.txt`: called cells that EmptyDrops dropped
- `filtered_species_calls.tsv`: every filtered barcode with its species call

In the one comparison run so far (lane 1, old symmetric 0.9 cutoff): all human and mixed calls were in the filtered set. 94% of pig calls were; the pig calls EmptyDrops dropped had low UMI counts (median 751). Only 15 filtered barcodes were left unassigned. **Re-run this on the current `_h0.8_p0.9_nd` output.**

---

## Step 5: Cross-species excess (`05_cross_species_excess.py`)

```bash
python3 05_cross_species_excess.py \
    --solo-dir      $D/Solo.out/GeneFull_Ex50pAS/raw \
    --exon-solo-dir $D/Solo.out/Gene/raw \
    --assignment    $D/classified_cells_h0.8_p0.9_nd/species_assignment.csv.gz \
    [--clusters pig_clusters.csv human_clusters.csv --barcode-filter '^lane1_'] \
    --outdir        $D/cross_species_excess_nd
```

This asks, in each direction (human reads in pig cells, pig reads in human cells), whether the other-species signal looks like **ambient RNA**, **mismapping**, or **something associated with the cell**.

### Why it tests pattern, not overall level
ρ is calculated *from* the other-species reads, so "total foreign UMIs vs expected" always matches by construction and tests nothing. Every test therefore looks at **which genes** and **which cells or clusters**.

### The tests

1. **Gene composition vs soup (does not use ρ).** Is gene *g* a bigger share of the foreign reads in cells than of that species' reads in the soup? This is a conditional binomial test on O_g given O_g + K_g, so noise in the soup counts is accounted for. `enrichment` = share in cells ÷ share in soup.
2. **Mismapping ceiling.** A host read that maps uniquely to a foreign gene (assembly gap, missing or diverged ortholog) shows up both in cells and in the soup, which is mostly host RNA. Its enrichment is always the same factor, whatever the mismapping rate:
   ```
   ceiling = (1 / foreign fraction in cells) × (foreign soup share / host soup share)
   ```
   That's about **6.2× for human genes in pig cells** and **about 117× for pig genes in human cells**. A mix of ambient RNA (1×) and mismapping falls between 1× and the ceiling and can't exceed it.

   | Enrichment | Interpretation |
   |---|---|
   | ≈ 1× | ambient RNA |
   | 1× to ceiling | ambiguous (partly mismapped, soup-vs-cell difference, or transfer) |
   | ≈ ceiling (`mismap_like`) | consistent with nearly all reads being mismapped host reads |
   | **significantly above the ceiling** (`q_above_ceiling` < 0.05) | not explained by ambient RNA or mismapping |
3. **Cells.** Foreign UMIs are modelled against cell size (negative binomial, log μ = a + b·log N). Per-cell q-values flag outlier cells. With `--exon-solo-dir`, excess cells whose foreign reads are also intronic are flagged `nuclear_like`; these are the same cells 03 calls `nuclear_doublet`.
4. **Clusters (with `--clusters`).**
   - **Observed vs expected foreign UMIs** per cluster, relative to the rest of the species.
   - **Cluster-vs-rest gene test** (`q_vs_rest`): artifacts present in every cluster cancel, so this is the column for transfer specific to a cell type.
   - **`q_nuclear_enriched`:** does the cluster have more nuclear-like cells than chance? (Fisher)
   - **`q_intronic`:** are the cluster's *ordinary* cells' foreign reads more intronic than the rest? (quasi-binomial)

   Cluster each species **on its own genes only**. If the other species' genes influence the clustering, finding a cluster "enriched for foreign RNA" becomes circular. Cluster barcodes are matched on the barcode sequence, so Seurat names like `lane1_<barcode>-1` work; use `--barcode-filter` to pick one lane from a merged object.

### Outputs
| File | Contents |
|---|---|
| `gene_excess.csv.gz` | One row per foreign gene × direction × group (`all` + each cluster): `O`, `K_soup`, `E`, `enrichment`, `p/q_up`, `p/q_down`, `mismap_ceiling`, `mismap_like`, `q_above_ceiling`; with clusters `O_rest`, `enrichment_vs_rest`, `p/q_vs_rest`; with exon counts `O_exonic`, `exonic_frac`, `soup_exonic_frac` |
| `cell_excess.csv.gz` | Per cell: `N`, `foreign`, `expected`, `ratio`, `p`, `q`; with exon counts `foreign_exonic_frac`, `intronic_q`, `nuclear_like`; `called_nuclear_doublet_03` |
| `cluster_excess.csv` | Per cluster: observed vs expected foreign UMIs, `q_up`, `rho_tracer`, exonic fraction, `q_intronic`, `n_nuclear_like`, `q_nuclear_enriched` |
| `summary.json` | Per direction: composition agreement, number of enriched / mismap-like / above-ceiling genes, size model, number of excess and nuclear-like cells, ambient UMIs per droplet |
| `gene_composition.png` | Observed vs soup-expected foreign UMIs per gene. Dashed line = 1×, dotted line = mismapping ceiling. Orange = mismap-like (hollow = above the line but not significantly), red = enriched below the ceiling, purple = significantly above the ceiling |
| `cell_foreign_vs_size.png` | Foreign UMIs vs cell size with the fitted trend. Red circles = excess, orange circles = excess + intronic |

`q` values are corrected within one direction × group, so compare them within a block. With huge counts, small effects become very significant; use `enrichment` to judge how big an effect is.

### Validation
- **Random cluster labels (null):** observed/expected = 0.99–1.01 for every cluster; 0 of about 62,000 cluster-vs-rest gene tests significant; no cluster intronic or nuclear-enriched.
- **Positive control:** labelling the nuclear doublets as a cluster gives `q_nuclear_enriched` ≈ 1e-293 and an exonic fraction of 0.42 vs 0.89.

---

## Results so far (pooled, no clusters yet)

| | Lane 1 | Lane 2 |
|---|---|---|
| Human / pig / mixed / nuclear-doublet cells | 1,140 / 18,364 / 456 / 113 | 1,137 / 18,068 / 443 / 117 |
| Human fraction of cells (input mix was 20% human) | 5.8% | 5.9% |
| Human share of soup | 7.1% | 7.1% |
| **ρ_pig / ρ_human** | **0.172 / 0.119** | **0.171 / 0.120** |
| Human genes enriched in pig cells vs soup / mismap-like / significantly above ceiling | 147 / 99 / **0** | 147 / 104 / **0** |
| Pig genes enriched in human cells / mismap-like / significantly above ceiling | 143 / 1 / **0** | 117 / 1 / **0** |
| Pig cells with excess human UMIs / of those, nuclear | 154 / 111 | 145 / 113 |
| Human exonic fraction in ordinary vs excess pig cells | 0.89 vs 0.54 | 0.90 vs 0.52 |
| Ambient per droplet, pig ÷ human droplets (same cell size) | 1.45× | 1.44× |

Interpretation:
- **Most human reads in pig cells are ambient RNA**: their gene composition matches the soup (Spearman 0.83). The most enriched genes (LINC00486, RPL6, EIF5, UBE2G1, SLC26A3, SOX2…) sit at the **mismapping ceiling**, so they're most likely pig reads mapping to human genes, not transfer. No gene is significantly above the ceiling in either direction.
- **Pig reads in human cells** show only 1–2× enrichment (ribosomal genes, GAPDH, HSP90AB1). That fits the usual difference between soup and cell-bound ambient RNA; no mismapping signature.
- **About 110 pig cells per lane (0.6%) carry human nuclear RNA**: excess human UMIs that are about half intronic, from long intron-rich genes (KIF26B, PTPRG, FBN2, WWOX, CADM1, MALAT1) and somite genes (FBN2, TCF12, EPHA4, GREB1L). The most likely explanation is a pig cell captured with a fragment of a human nucleus. **Pig–human cell fusion** is the less likely alternative: the excess is only about 70 UMIs per cell, far less than a working nucleus. These cells are the best candidates for follow-up and are kept out of the clean pig set.
- Pig droplets carry about 1.45× more ambient RNA than human droplets at the same size. That's a hint, not evidence of transfer: the two cell types may differ in fragility or surface-bound debris.
- About 25 pig cells near the `--min-umi` cutoff have excess human UMIs that are **not** intronic. They're probably low-quality droplets.

---

## Step 6: Per-species analysis (single-genome alignment → SoupX → subset → Seurat)

Gene counts for each species come from aligning the **whole library** to that species' genome alone (`run_aligner.sh` with the human or pig index): `shared_mount/ChimeraTLS_20H80P_mapped2humanonly_lane{1,2}` and `…mapped2pigonly_lane{1,2}`. Each species gets its own gene models (e.g. the 3′-extended pig GTF), and reads from conserved genes aren't lost as human/pig multimappers.

`run_seuratv5.sh` then runs SoupX with the ρ from step 3 and subsets to that species' barcodes. Add two keys to the `fastq_alignment` sample in `config.json` (one species per run):

```json
"SAMPLE_1": {
    "NAME": "ChimeraTLS_20H80P_mapped2pigonly",
    "species": "pig",
    "COMBINE_LANES": true,
    "SPECIES_BARCODES": "/w5home/bmoore/scRNAseq/LiFangChu/10x_somite_data/chimera_data/ChimeraTLS_20H80P_3pext/<03 outdir>/barcodes_pig.txt",
    "SOUPX_RHO": 0.17,
    ...
}
```

- **`species`** must be the species kept (`pig` or `human`), not `human-pig`. Downstream steps use it for mitochondrial genes and cell-cycle genes.
- **`SPECIES_BARCODES`** is `barcodes_pig.txt` or `barcodes_human.txt` from step 3, as an absolute path. Mixed and nuclear-doublet barcodes are already left out of these files.
- **`SOUPX_RHO`** is `rho_pig` or `rho_human` from step 3's `summary.json`. Use the ρ for the species kept.
- With `COMBINE_LANES: true`, both keys are **single values** (the launcher rejects per-lane maps). For per-lane folders (`COMBINE_LANES: false`) they can be per-lane maps such as `{"lane1": …, "lane2": …}`, matched by the suffix after `NAME_`. Don't merge per-lane objects of the same library in the Seurat step, though: see Step 2.
- Use the barcode files and ρ from the **same** step-3 run, and update both when you re-run step 3.

Run with `source run_seuratv5.sh config.json`. The launcher checks these keys before starting: files exist, ρ is between 0 and 1, the species is valid, and lanes match. It mounts each barcode folder read-only into the container.

What happens per lane (`sc_pipeline/src/process_lane.R`):
1. **SoupX runs on the whole library.** The soup profile comes from all empty droplets in that alignment, so other-species ambient RNA mapped onto this genome is included. The cell set is STARsolo's filtered cells **plus** any species barcodes EmptyDrops missed: in the comparison run so far, about 6% of pig calls weren't in STARsolo's filtered set.
2. **`setContaminationFraction(rho)` replaces `autoEstCont()`.** SoupX is still needed even though ρ is known: ρ says *how much* to remove, and `adjustCounts()` removes it gene by gene, including the same-species ambient RNA that the tracer can't see.
3. **The corrected counts are subset to `SPECIES_BARCODES`.** Everything downstream (filtering, scDblFinder, normalization, clustering, annotation) sees only that species' cells.
4. A per-lane summary goes to `output/output_<timestamp>/species_subset_<lane>.tsv`: barcodes given, how many are missing from the raw matrix, how many were added beyond STARsolo's filtered cells, cells kept, and median UMIs before and after SoupX.

Without these keys the pipeline behaves as before (STARsolo filtered cells, `autoEstCont()`).

## Next steps

1. **Combine the lanes** (Step 2) for all three alignments: combined reference, human-only and pig-only. Then re-run 03 with one `--min-umi`, plus 04 and 05, on the combined reference, and update `SPECIES_BARCODES` and `SOUPX_RHO` to single values with `COMBINE_LANES: true`.
2. **Run step 6 for each species** (pig and human, as separate runs).
3. **Cluster each species on its own genes**, then re-run step 5 with `--clusters` to test per-cluster and cell-type-specific transfer.
4. **For any candidate gene or cell**, check reads in the BAM: mismatches (`nM`), whether the read also aligns almost as well to the other genome, and whether the host assembly has the gene.
5. **Validation outside sequencing:** tdTomato protein in pig cells (flow cytometry or imaging), and RNA FISH for a human-specific transcript.

---

## Known issues and gotchas

- **Re-run 03 on any output made before 2026-11-18.** An earlier version of 03 read the matrix size line as a count. That created a fake "pig cell" (the last barcode, with about 60M UMIs), which appeared in `barcodes_pig.txt` and pulled ρ_pig down to about 0.10. 04 refuses to run on such output.
- **Rebuilding the reference:** `01_build_combined_reference.sh` reproduces the reference on disk (rebuilt 2026-11-12) from the local source files. It refuses to overwrite an existing index unless `FORCE=1`. Before indexing it checks that every GTF sequence name exists in the FASTA, and after indexing it checks STAR's log and the per-species gene counts.
- **Use one feature matrix throughout** (`GeneFull_Ex50pAS`) for 03, 04 and 05. Numbers shift slightly between `GeneFull` and `GeneFull_Ex50pAS`.
- **Species fraction:** only about 6% of cells are human, although the culture was mixed at 20% human. That could be differential growth or survival by day 3, or loss of human cells during dissociation. Human cells also carry more ambient RNA than pig cells, which may point to fragility.
- **Mismapped genes in the human annotation** (LINC00486, RPL6, EIF5, SOX2, SOX2-OT…) will also show low-level "expression" in pig cells in any combined-reference analysis. Don't read these as human expression.
