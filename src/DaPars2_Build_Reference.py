#!/usr/bin/env python3
"""Build a DaPars2 reference from selected genes and RNA-seq-supported terminal exons.

Created: 2026-10-23
Updated: 2026-10-23
Inputs: GTF, gene-type and gene-ID lists, and a raw BigWig/mapped-alignment sample TSV.
Outputs: DaPars2 BED, standard-coordinate IGV BED, model/region metadata and summary.
Workflow: Gene selection -> full terminal-exon coverage -> same-5'-boundary collapse.
Key parameters: Coverage >=2 per 50M alignments; >=75% bases; >=20% samples.
Assumptions: Unstranded raw coverage and annotation use the same genome assembly.
QC: Preserve gene/transcript versions, gene/locus identity and exact covered-base counts.
Run: python src/DaPars2_Build_Reference.py --help
"""

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
from fractions import Fraction
import gzip
import math
from pathlib import Path
import re

import numpy as np
import pyBigWig


ATTRIBUTES = re.compile(r'(\w+)\s+"([^"]*)"')
ENSEMBL_VERSION = re.compile(r'^(ENS[A-Za-z0-9]*G\d+)\.\d+(?:_[A-Za-z0-9]+)*$')
WORKER_INTERVALS = None
WORKER_MIN_COVERAGE = None
WORKER_NORMALIZATION_DEPTH = None


def stable_gene_id(gene_id):
    """Match versioned GENCODE/Ensembl gene lists without changing output IDs."""
    match = ENSEMBL_VERSION.fullmatch(gene_id)
    return match.group(1) if match else gene_id


def chromosome_key(chrom):
    return tuple(int(part) if part.isdigit() else part for part in re.split(r'(\d+)', chrom))


def read_list(path):
    values = set()
    with path.open() as source:
        for number, line in enumerate(source, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if len(line.split()) != 1:
                raise ValueError(f'{path}:{number}: expected one value per line, without a header')
            values.add(line)
    if not values:
        raise ValueError(f'{path}: selection list is empty')
    return values


def read_samples(path):
    samples, identifiers, paths = [], set(), set()
    with path.open() as source:
        reader = csv.DictReader(source, delimiter='\t')
        if not {'sample_id', 'bigwig', 'mapped_alignments'} <= set(reader.fieldnames or []):
            raise ValueError('Sample TSV requires sample_id, bigwig and mapped_alignments columns')
        for row in reader:
            if not any(row.values()):
                continue
            sample_id = row['sample_id'].strip()
            bigwig = Path(row['bigwig'].strip()).expanduser()
            if not bigwig.is_absolute():
                bigwig = path.parent / bigwig
            bigwig = bigwig.resolve()
            mapped = float(row['mapped_alignments'])
            if not math.isfinite(mapped) or mapped <= 0 or not mapped.is_integer():
                raise ValueError(f'{sample_id}: mapped_alignments must be a positive integer')
            # Duplicate libraries would inflate the cohort support fraction silently.
            if not sample_id or sample_id in identifiers or bigwig in paths:
                raise ValueError(f'{sample_id}: sample identifiers and BigWig paths must be unique')
            samples.append(dict(sample_id=sample_id, bigwig=str(bigwig), mapped_alignments=int(mapped)))
            identifiers.add(sample_id)
            paths.add(bigwig)
    if not samples:
        raise ValueError('Sample TSV has no samples')
    return samples


def read_annotation(path, gene_types, gene_ids, chromosomes=None):
    """Select genes and their complete last exon in transcriptional order."""
    genes, transcripts = {}, {}
    opener = gzip.open if path.suffix in {'.gz', '.bgz'} else open
    with opener(path, 'rt') as source:
        for number, line in enumerate(source, 1):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{number}: expected nine GTF columns')
            chrom, feature, strand = fields[0], fields[2], fields[6]
            if feature not in {'gene', 'exon'} or (chromosomes and chrom not in chromosomes):
                continue
            attrs = dict(ATTRIBUTES.findall(fields[8]))
            gene_id = attrs['gene_id']
            key = (gene_id, chrom, strand)
            gene = genes.setdefault(key, dict(gene_id=gene_id, gene_symbol='', gene_type=''))
            for name, value in (('gene_symbol', attrs.get('gene_name', '')),
                                ('gene_type', attrs.get('gene_type') or attrs.get('gene_biotype', ''))):
                if value and (not gene[name] or feature == 'gene'):
                    gene[name] = value
            if feature != 'exon' or stable_gene_id(gene_id) not in gene_ids:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            if strand not in {'+', '-'} or start < 0 or end <= start:
                raise ValueError(f'{path}:{number}: selected exons require +/-, start>=1 and end>=start')
            transcript_id = attrs['transcript_id']
            transcript_key = key + (transcript_id,)
            previous = transcripts.get(transcript_key)
            if previous is None or (strand == '+' and start > previous[0]) or (strand == '-' and start < previous[0]):
                transcripts[transcript_key] = (start, end)

    missing_types = [key for key, gene in genes.items()
                     if stable_gene_id(key[0]) in gene_ids and not gene['gene_type']]
    if missing_types:
        raise ValueError(f'{missing_types[0][0]}: gene_type/gene_biotype is absent from the GTF')
    typed = {key for key, gene in genes.items() if gene['gene_type'] in gene_types}
    selected = {key for key in typed if stable_gene_id(key[0]) in gene_ids}
    exact_models = defaultdict(set)
    for key, interval in transcripts.items():
        if key[:3] in selected:
            exact_models[key[:3] + interval].add(key[3])
    models = []
    for (gene_id, chrom, strand, start, end), transcript_ids in exact_models.items():
        gene = genes[(gene_id, chrom, strand)]
        symbol = gene['gene_symbol'] or gene_id
        if any('|' in value for value in [gene_id, symbol, *transcript_ids]) or any(',' in value for value in transcript_ids):
            raise ValueError(f'{gene_id}: identifiers/symbols cannot contain |; transcript IDs cannot contain commas')
        models.append(dict(gene_id=gene_id, stable_gene_id=stable_gene_id(gene_id), gene_symbol=symbol,
                           gene_type=gene['gene_type'], chrom=chrom, strand=strand, start0=start, end0=end,
                           anchor0=start if strand == '+' else end, length=end-start,
                           transcript_ids=','.join(sorted(transcript_ids)), n_transcripts=len(transcript_ids)))
    models.sort(key=lambda row: (chromosome_key(row['chrom']), row['start0'], row['end0'], row['gene_id'], row['strand']))
    if not models:
        raise ValueError('No terminal-exon models selected; check gene IDs, gene types and chromosomes')
    stats = dict(n_annotation_gene_loci=len(genes), n_gene_loci_after_type_selection=len(typed),
                 n_gene_loci_after_gene_list=len(selected),
                 n_requested_gene_ids=len(gene_ids),
                 n_unmatched_gene_ids=len(gene_ids - {stable_gene_id(key[0]) for key in genes}),
                 n_selected_transcripts=sum(row['n_transcripts'] for row in models),
                 n_terminal_models=len(models))
    return models, stats


def required_count(size, fraction):
    """Compute ceil(size*fraction) exactly, including decimal boundaries such as 7%."""
    return math.ceil(size * Fraction(str(fraction)))


def prepare_intervals(models, base_fraction):
    geometries = sorted({(row['chrom'], row['start0'], row['end0']) for row in models},
                        key=lambda row: (chromosome_key(row[0]), row[1], row[2]))
    indexes = {geometry: i for i, geometry in enumerate(geometries)}
    for row in models:
        row['geometry_index'] = indexes[(row['chrom'], row['start0'], row['end0'])]
    by_chrom = defaultdict(list)
    for index, (chrom, start, end) in enumerate(geometries):
        by_chrom[chrom].append((index, start, end, required_count(end-start, base_fraction)))
    return {chrom: np.asarray(rows, dtype=np.int64) for chrom, rows in by_chrom.items()}, len(geometries)


def covered_to(points, intervals):
    """Count qualifying bases before each boundary without allocating per-base coverage."""
    if len(intervals) == 0:
        return np.zeros(len(points), dtype=np.int64)
    widths = intervals[:, 1] - intervals[:, 0]
    cumulative = np.r_[0, np.cumsum(widths)]
    index = np.searchsorted(intervals[:, 0], points, side='right') - 1
    safe = np.maximum(index, 0)
    partial = np.clip(points - intervals[safe, 0], 0, widths[safe])
    return np.where(index < 0, 0, cumulative[safe] + partial)


def initialize_worker(intervals, min_coverage, normalization_depth):
    global WORKER_INTERVALS, WORKER_MIN_COVERAGE, WORKER_NORMALIZATION_DEPTH
    WORKER_INTERVALS = intervals
    WORKER_MIN_COVERAGE = min_coverage
    WORKER_NORMALIZATION_DEPTH = normalization_depth


def measure_sample(sample):
    size = sum(len(rows) for rows in WORKER_INTERVALS.values())
    passing = np.zeros(size, dtype=np.int32)
    raw_cutoff = WORKER_MIN_COVERAGE * sample['mapped_alignments'] / WORKER_NORMALIZATION_DEPTH
    with pyBigWig.open(sample['bigwig']) as bigwig:
        lengths = bigwig.chroms()
        for chrom, rows in WORKER_INTERVALS.items():
            # Missing chromosomes cannot be distinguished from an assembly/naming mismatch.
            if chrom not in lengths or int(rows[:, 2].max()) > lengths[chrom]:
                raise ValueError(f"{sample['sample_id']}: annotation coordinates exceed or are absent from BigWig chromosome {chrom}")
            runs = bigwig.intervals(chrom, int(rows[:, 1].min()), int(rows[:, 2].max()))
            if runs is None:
                continue
            runs = np.asarray(runs, dtype=np.float64)
            if np.any(~np.isfinite(runs[:, 2])) or np.any(runs[:, 2] < 0):
                raise ValueError(f"{sample['sample_id']}: expected finite, nonnegative raw read coverage")
            good = runs[runs[:, 2] >= raw_cutoff, :2].astype(np.int64)
            counts = covered_to(rows[:, 2], good) - covered_to(rows[:, 1], good)
            passing[rows[:, 0]] = counts >= rows[:, 3]
    return passing


def measure_support(samples, intervals, n_geometries, min_coverage, normalization_depth, workers):
    support = np.zeros(n_geometries, dtype=np.int32)
    workers = min(workers, len(samples))
    if workers == 1:
        initialize_worker(intervals, min_coverage, normalization_depth)
        for sample in samples:
            support += measure_sample(sample)
    else:
        with ProcessPoolExecutor(max_workers=workers, initializer=initialize_worker,
                                 initargs=(intervals, min_coverage, normalization_depth)) as pool:
            for processed, passing in enumerate(pool.map(measure_sample, samples), 1):
                support += passing
                if processed % 25 == 0 or processed == len(samples):
                    print(f'Coverage: {processed}/{len(samples)} samples', flush=True)
    return support


def collapse_models(models):
    """Merge only passing models of the same gene/locus and exact transcriptional 5' anchor."""
    groups = defaultdict(list)
    for row in models:
        if row['coverage_pass']:
            groups[(row['gene_id'], row['chrom'], row['strand'], row['anchor0'])].append(row)
    regions, excluded = [], 0
    for group in groups.values():
        longest = max(group, key=lambda row: row['length'])
        transcripts = sorted({value for row in group for value in row['transcript_ids'].split(',')})
        region = dict(longest, transcript_ids=','.join(transcripts), n_transcripts=len(transcripts),
                      n_passing_models=len(group), dapars2_start=longest['start0']+1, dapars2_end=longest['end0'])
        # A one-base exon cannot be represented in the upstream start0+1,end0 convention.
        if region['dapars2_end'] <= region['dapars2_start']:
            excluded += 1
            continue
        region['name'] = '|'.join([region['gene_id'], region['gene_symbol'], region['chrom'], region['strand'], region['transcript_ids']])
        regions.append(region)
    regions.sort(key=lambda row: (chromosome_key(row['chrom']), row['start0'], row['end0'], row['name']))
    return regions, excluded


def write_table(path, rows, columns):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt') as output:
        writer = csv.DictWriter(output, columns, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def build_reference(args):
    types = read_list(args.gene_types)
    gene_ids = {stable_gene_id(value) for value in read_list(args.gene_list)}
    samples = read_samples(args.samples)
    models, stats = read_annotation(args.gtf, types, gene_ids, set(args.chromosomes) if args.chromosomes else None)
    intervals, n_geometries = prepare_intervals(models, args.min_base_fraction)
    print(f'Selected {stats["n_gene_loci_after_gene_list"]} gene loci; {len(models)} terminal-exon models; '
          f'{len(samples)} samples; workers={min(args.workers, len(samples))}', flush=True)
    support = measure_support(samples, intervals, n_geometries, args.min_coverage, args.normalization_depth, args.workers)
    minimum_samples = required_count(len(samples), args.min_sample_fraction)
    for row in models:
        row.update(n_required_bases=required_count(row['length'], args.min_base_fraction), n_samples=len(samples),
                   minimum_samples=minimum_samples, n_pass=int(support[row['geometry_index']]))
        row['pass_fraction'] = row['n_pass'] / len(samples)
        row['coverage_pass'] = row['n_pass'] >= minimum_samples
    regions, excluded = collapse_models(models)
    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    model_columns = ['gene_id', 'stable_gene_id', 'gene_symbol', 'gene_type', 'chrom', 'strand', 'start0', 'end0',
                     'anchor0', 'length', 'transcript_ids', 'n_transcripts', 'n_required_bases', 'n_samples',
                     'minimum_samples', 'n_pass', 'pass_fraction', 'coverage_pass']
    write_table(Path(str(prefix)+'.terminal_models.tsv.gz'), models, model_columns)
    write_table(Path(str(prefix)+'.regions.tsv.gz'), regions,
                model_columns + ['n_passing_models', 'dapars2_start', 'dapars2_end', 'name'])
    for suffix, start, end in [('.dapars2.bed', 'dapars2_start', 'dapars2_end'),
                               ('.full_exon.igv.bed', 'start0', 'end0')]:
        with Path(str(prefix)+suffix).open('w') as output:
            writer = csv.writer(output, delimiter='\t', lineterminator='\n')
            writer.writerows((row['chrom'], row[start], row[end], row['name'], 0, row['strand']) for row in regions)
    stats.update(n_samples=len(samples), n_unique_geometries=n_geometries, minimum_samples=minimum_samples,
                 min_coverage=args.min_coverage, normalization_depth=args.normalization_depth,
                 min_base_fraction=args.min_base_fraction, min_sample_fraction=args.min_sample_fraction,
                 n_models_passing_coverage=sum(row['coverage_pass'] for row in models),
                 n_regions_excluded_coordinate_convention=excluded, n_final_regions=len(regions),
                 n_final_gene_loci=len({(row['gene_id'],row['chrom'],row['strand']) for row in regions}))
    write_table(Path(str(prefix)+'.summary.tsv'), [dict(metric=key, value=value) for key, value in stats.items()], ['metric', 'value'])
    print(f'Finished: {len(regions)} reference regions; {excluded} one-base regions excluded by DaPars2 coordinates.', flush=True)
    return stats


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--gtf', type=Path, required=True, help='GTF annotation (.gtf, .gtf.gz or .gtf.bgz)')
    parser.add_argument('--gene-types', type=Path, required=True, help='One retained gene_type/gene_biotype per line')
    parser.add_argument('--gene-list', type=Path, required=True, help='One externally selected gene_id per line')
    parser.add_argument('--samples', type=Path, required=True, help='TSV: sample_id, bigwig, mapped_alignments')
    parser.add_argument('--output-prefix', type=Path, required=True, help='Output path prefix')
    parser.add_argument('--chromosomes', nargs='+', help='Optional exact GTF/BigWig sequence names; default: all selected genes')
    parser.add_argument('--min-coverage', type=float, default=2, help='Minimum normalized per-base coverage (default: 2)')
    parser.add_argument('--normalization-depth', type=float, default=50_000_000, help='Reference mapped alignments (default: 50000000)')
    parser.add_argument('--min-base-fraction', type=float, default=0.75, help='Minimum fraction of full-exon bases meeting coverage (default: 0.75)')
    parser.add_argument('--min-sample-fraction', type=float, default=0.20, help='Minimum cohort fraction supporting a model (default: 0.20)')
    parser.add_argument('--workers', type=int, default=4, help='Sample workers (default: 4)')
    args = parser.parse_args()
    if args.workers < 1 or not all(math.isfinite(value) and value > 0 for value in (args.min_coverage, args.normalization_depth)):
        parser.error('workers, min-coverage and normalization-depth must be positive and finite')
    if not all(math.isfinite(value) and 0 < value <= 1 for value in (args.min_base_fraction, args.min_sample_fraction)):
        parser.error('base/sample fractions must be in (0,1]')
    try:
        build_reference(args)
    except (ValueError, OSError, RuntimeError, KeyError) as error:
        parser.exit(2, f'error: {error}\n')


if __name__ == '__main__':
    main()
