# RNA-seq-guided terminal-exon reference preparation

The new entry point is `src/DaPars2_Build_Reference.py`. It complements `src/DaPars_Extract_Anno.py`; it does not change DaPars2 fitting or PDUI calculation.

## Inputs

Supply an uncompressed or gzip-compressed GTF in the same assembly and with the same chromosome names as the BigWigs. Exons need `gene_id` and `transcript_id`. Gene type is read from `gene_type` or `gene_biotype`, on gene or exon records. `gene_name` supplies the symbol; if absent, the original gene ID is used. Explicit gene records are optional.

All transcript biotypes within a selected gene are eligible. Terminal exons are identified in transcriptional order: greatest exon start on `+`, smallest exon start on `-`. These are complete exons and can include coding sequence; the tool does not restrict them to annotated 3′UTR features or a fixed-length window.

Chromosomes are not hard-coded to human autosomes. By default all selected gene loci are considered; use, for example, `--chromosomes chr1 chr2 chr21` to restrict them. Coordinate mismatches or absent BigWig chromosomes are errors rather than assumed zero coverage. Missing intervals on a chromosome that is present count as zero.

### Selection lists

Both files are required, one value per line, without a header. Blank lines and lines beginning with `#` are ignored.

`gene_types.txt` specifies retained annotation types:

```text
protein_coding
```

`gene_list.txt` specifies retained gene IDs:

```text
ENSG00000000419
ENSG00000001036
```

The gene list is an external inclusion list. Generate it using an expression criterion appropriate to your study; no TPM/count matrix or eQTL-specific expression rule is built in. Both gene-type and gene-ID selections must pass.

GENCODE/Ensembl gene versions are ignored for list matching, including lift37 suffixes: `ENSG00000000419`, `ENSG00000000419.14`, and `ENSG00000000419.14_13` match the same stable identifier. Other namespaces are matched exactly. Outputs retain the original GTF gene and transcript IDs. Gene identity, chromosome and strand remain separate even when coordinates overlap or an identifier occurs at multiple loci.

### Samples and normalization

The sample file is a tab-separated table with a header:

```text
sample_id	bigwig	mapped_alignments
sample_1	coverage/sample_1.bw	50000000
sample_2	coverage/sample_2.bw	100000000
```

Relative BigWig paths are relative to the sample table's directory. Each sample ID and BigWig path must be unique. One invocation treats all rows as one cohort; separate tissues/cohorts can use separate lists, sample tables and output prefixes.

BigWigs must contain raw, unnormalized, nonnegative, unstranded read coverage. They must not be CPM/RPKM/binned/smoothed or signed strand tracks. Overlapping mates count separately. `mapped_alignments` counts alignment records from the same source BAM/library, rather than read pairs or fragments; use a consistent definition across samples. To reproduce the original project normalization, supply the whole-genome mapped count even when evaluating autosomal regions.

The default normalization is:

```text
normalized_depth = raw_depth * 50000000 / mapped_alignments
```

The raw cutoff is therefore `2 * mapped_alignments / 50000000`: a 50M-alignment sample needs raw coverage ≥2, and a 100M-alignment sample needs ≥4. Fractional signals are not rounded before comparison. This reference-selection normalization is independent of the depth normalization performed later by DaPars2.

Starting from an already filtered, coordinate-sorted BAM, example raw-track preparation is:

```bash
samtools view -c -F 4 sample.bam
bedtools genomecov -bga -split -ibam sample.bam > sample.bedGraph
LC_ALL=C sort -k1,1 -k2,2n sample.bedGraph > sample.sorted.bedGraph
bedGraphToBigWig sample.sorted.bedGraph chrom.sizes sample.bw
```

Use the mapped count in the sample table. The BAM, chromosome sizes and GTF must describe the same assembly. The reference builder itself only requires NumPy and pyBigWig; these external tools prepare coverage inputs.

## Selection and collapse

1. Select genes using the type and gene-ID files, then identify each transcript's complete terminal exon.
2. Share coverage measurement for identical genomic geometries without merging their gene/strand identities.
3. Per sample, count bases with normalized coverage ≥`--min-coverage`. Require at least `ceil(exon_length * min_base_fraction)` qualifying bases.
4. Retain a model if at least `ceil(sample_count * min_sample_fraction)` samples pass. Comparisons are inclusive; at 20%, both 202 and 204 samples need 41 supporters, while six samples need two.
5. Group only passing models by gene, chromosome, strand and exact transcriptional 5′ boundary: `start0` on `+`, `end0` on `-`. Keep the longest passing model and the union of passing transcript IDs. Failed models do not contribute IDs.

The default thresholds are a study recipe, not universal evidence of transcript use. Coverage supports annotated models but does not establish a true cleavage site or uniquely assign overlapping, unstranded reads to a gene. Expressed-gene selection can reduce interference but cannot remove it when both overlapping genes are expressed.

Unlike the legacy extractor, this entry point does not run opposite-strand subtraction or otherwise trim a passing exon. DaPars2's later trimming, minimum-length and fitting rules remain unchanged; passing the reference filter does not guarantee a fitted APA event.

## Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--min-coverage` | `2` | Normalized per-base coverage cutoff |
| `--normalization-depth` | `50000000` | Mapped alignments for reference normalization |
| `--min-base-fraction` | `0.75` | Required fraction of the complete exon |
| `--min-sample-fraction` | `0.20` | Required fraction of supplied samples |
| `--workers` | `4` | Sample processes, bounded by sample count |
| `--chromosomes` | all selected loci | Optional exact sequence names, separated by spaces |

Choose workers according to available CPU and memory. Coverage is integrated over BigWig runs instead of allocating chromosome-sized per-base arrays; each worker holds one sample at a time.

## Outputs and coordinates

For `--output-prefix results/reference`, five files are produced:

- `reference.dapars2.bed`: BED6-like input using the upstream extractor's **`start0+1,end0`** convention. Point `Annotated_3UTR` here.
- `reference.full_exon.igv.bed`: identical regions as standard **0-based, half-open `[start0,end0)`** BED6 for IGV/bedtools.
- `reference.regions.tsv.gz`: final coordinates, name, gene type, original IDs, transcript membership, passing-model count and support of the longest passing model.
- `reference.terminal_models.tsv.gz`: all selected exact exon models, including failed coverage decisions; includes required bases, sample counts, supporters and retention flags.
- `reference.summary.tsv`: selection counts, unmatched gene IDs within selected chromosomes, thresholds and excluded one-base regions.

A one-base region cannot have positive span in the legacy DaPars2 convention; it is excluded after collapse and counted. No other reference-length threshold is added.

The ID is `gene_id|gene_symbol|chromosome|strand|transcript_ids`. Full gene/transcript versions are retained; merged transcripts are comma-separated. Input IDs/symbols cannot contain `|`, and individual transcript IDs cannot contain commas, because these characters delimit the output name. Quote names when passing them through a shell, because `|` is a shell operator. Coordinates and support remain separate metadata, not extra ID fields.
