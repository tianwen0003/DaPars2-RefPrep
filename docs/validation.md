# Validation

Validated on 2026-10-23 against the original project implementation and the unchanged upstream DaPars2 quantifier. These checks establish implementation agreement within the tested scope; they do not establish the biological accuracy of an annotation or APA call.

## Reproducible synthetic checks

Run from the repository root after installing `requirements-reference.txt`:

```bash
python -m unittest discover -s tests -v
```

All **10 tests passed** in both tested environments:

| Python | NumPy | pyBigWig |
| --- | --- | --- |
| 3.9.7 | 1.23.5 | 0.3.18 |
| 3.11.11 | 1.26.4 | 0.3.23 |

The tests create explicitly synthetic GTFs and BigWigs in temporary directories. They cover full terminal-exon selection on both strands, gene-type and external gene-list selection, version preservation, separate genes with identical coordinates, exact 5′-anchor collapse, exclusion of failed transcript models, normalization by mapped alignments, inclusive base and sample thresholds, integer ceiling boundaries, both output coordinate conventions, and serial/parallel agreement. They also exercise gene metadata supplied after exon records, `gene_biotype` annotations without gene records, arbitrary dotted IDs, duplicate libraries, missing chromosomes, negative signal, and reserved transcript delimiters.

## Real-cohort comparison

The builder was run separately on **chr21** for two cohorts using GENCODE v48lift37 basic, externally selected expressed protein-coding genes, raw BigWigs, and the default 2/75%/20% rule. Each cohort used two workers. Results were compared with the chr21 portion of the existing formal project references.

| Cohort | Samples | Terminal-exon models | Passing models | Collapsed regions | Final gene loci | Required supporting samples |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Cartilage | 204 | 393 | 306 | 171 | 121 | 41 |
| Synovium | 202 | 379 | 320 | 175 | 121 | 41 |

All selected model coordinates, transcript memberships, required base counts, per-model sample-support counts and retention decisions matched. Collapsed coordinates, complete IDs, passing-model memberships and support metadata also matched. Both the native DaPars2 BED and standard-coordinate IGV BED were **textually identical** to their corresponding formal chr21 reference rows.

The private cohort inputs and generated test outputs are not distributed in this repository. Reproduce this comparison with the relevant external selection lists and sample tables using `--chromosomes chr21`. This was a chromosome-level check, not a whole-genome or multi-annotation benchmark.

## Upstream quantification compatibility

An additional smoke test ran the unchanged `src/DaPars2_Multi_Sample_Multi_Chr.py` CLI on three regions from the new chr21 reference. The input covered both strands, two distinct regions from the same gene, complete versioned IDs and merged transcript lists. Coverage was explicitly **synthetic uniform depth 40**, with one sample and one worker; it was not the real-cohort RNA-seq signal.

The quantifier exited successfully, wrote all three fitted rows, and preserved their full names, using Python 3.11.11, NumPy 1.26.4 and SciPy 1.14.1. This checks the BED/identifier interface and execution path. Uniform synthetic coverage does not validate proximal-site accuracy or PDUI accuracy. The quantifier's own coverage, trimming and fitting rules remain unchanged.

The original three source scripts, test-dataset archive and license are unchanged from the fork's starting commit, `fb81c6c`. Its original README is preserved in [UPSTREAM_README.md](UPSTREAM_README.md).
