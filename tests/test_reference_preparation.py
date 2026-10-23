"""Synthetic, deterministic checks of reference selection, coordinates and cohort support."""

import csv
import gzip
import importlib.util
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pyBigWig


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'src/DaPars2_Build_Reference.py'
SPEC = importlib.util.spec_from_file_location('refprep', SCRIPT)
refprep = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(refprep)


class ReferencePreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='dapars2-refprep-test-')
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.gtf = self.directory / 'annotation.gtf.gz'
        self.types = self.directory / 'gene_types.txt'
        self.genes = self.directory / 'gene_list.txt'
        self.samples = self.directory / 'samples.tsv'
        self.types.write_text('# retained gene types\nprotein_coding\n')
        self.genes.write_text('\n'.join(f'ENSG{i:011d}' for i in (1, 2, 3, 4, 6, 7))+'\n')
        transcripts = [
            (1, '+', 'A1.1', [(0, 20), (100, 200)]),
            (1, '+', 'A1_DUP.4', [(100, 200)]),
            (1, '+', 'A2.2', [(0, 20), (100, 240)]),
            (1, '+', 'A3.3', [(120, 220)]),
            (1, '+', 'A_BAD.1', [(100, 300)]),
            (2, '-', 'B1.1', [(400, 500), (700, 720)]),
            (2, '-', 'B2.2', [(360, 500), (700, 720)]),
            (2, '-', 'B_BAD.1', [(300, 500)]),
            (3, '+', 'C1.1', [(100, 240)]),
            (4, '+', 'NONCODING.1', [(100, 240)]),
            (5, '+', 'NOT_LISTED.1', [(100, 240)]),
            (6, '-', 'ONE_BASE.1', [(20, 21)]),
            (7, '+', 'G100.1', [(800, 900)]),
            (7, '+', 'G101.1', [(800, 901)]),
        ]
        with gzip.open(self.gtf, 'wt') as output:
            for number, strand, transcript, exons in transcripts:
                chrom = 'chr2' if number == 6 else 'chr1'
                gene_type = 'lncRNA' if number == 4 else 'protein_coding'
                attrs = f'gene_id "ENSG{number:011d}.2_7"; gene_name "GENE{number}"; gene_type "{gene_type}"; transcript_id "{transcript}";'
                for start, end in exons:
                    output.write(f'{chrom}\ttest\texon\t{start+1}\t{end}\t.\t{strand}\t.\t{attrs}\n')
        for number in range(6):
            value = [2, 4, 1, 0, 0, 0][number]
            self.write_bigwig(self.directory / f's{number}.bw', value)
        self.samples.write_text('sample_id\tbigwig\tmapped_alignments\n'+''.join(
            f's{i}\ts{i}.bw\t{100000000 if i == 1 else 50000000}\n' for i in range(6)))

    def write_bigwig(self, path, value, chromosomes=('chr1', 'chr2')):
        with pyBigWig.open(str(path), 'w') as output:
            output.addHeader([(chrom, 1000) for chrom in chromosomes])
            if value and 'chr1' in chromosomes:
                output.addEntries(['chr1']*3, [100, 360, 800], ends=[240, 500, 875], values=[float(value)]*3)
            if value and 'chr2' in chromosomes:
                output.addEntries(['chr2'], [20], ends=[21], values=[float(value)])

    def run_cli(self, prefix='reference', workers=1, success=True, extra=()):
        command = [sys.executable, str(SCRIPT), '--gtf', str(self.gtf), '--gene-types', str(self.types),
                   '--gene-list', str(self.genes), '--samples', str(self.samples),
                   '--output-prefix', str(self.directory/prefix), '--workers', str(workers), *extra]
        result = subprocess.run(command, capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def table(self, suffix, prefix='reference'):
        path = self.directory / (prefix+suffix)
        opener = gzip.open if path.suffix == '.gz' else open
        with opener(path, 'rt') as source:
            return list(csv.DictReader(source, delimiter='\t'))

    def test_full_exons_identity_versions_and_passing_transcripts(self):
        self.run_cli()
        models = self.table('.terminal_models.tsv.gz')
        regions = self.table('.regions.tsv.gz')
        self.assertEqual(len(models), 11)
        self.assertEqual(len(regions), 5)
        self.assertTrue(all(row['n_pass'] == '2' for row in regions))
        self.assertTrue(all(row['minimum_samples'] == '2' for row in regions))
        self.assertTrue(all('.2_7|' in row['name'] for row in regions))
        alpha = [row for row in regions if row['gene_symbol'] == 'GENE1']
        self.assertEqual({(row['start0'], row['end0']) for row in alpha}, {('100', '240'), ('120', '220')})
        merged = next(row for row in alpha if row['start0'] == '100')
        self.assertEqual(merged['transcript_ids'], 'A1.1,A1_DUP.4,A2.2')
        self.assertEqual(merged['n_passing_models'], '2')
        beta = next(row for row in regions if row['gene_symbol'] == 'GENE2')
        self.assertEqual((beta['start0'], beta['end0'], beta['anchor0']), ('360', '500', '500'))
        # Identical geometry belonging to another gene remains another phenotype.
        same_geometry = [row for row in regions if (row['start0'], row['end0']) == ('100', '240')]
        self.assertEqual(len(same_geometry), 2)
        bad = [row for row in models if 'BAD' in row['transcript_ids'] or row['transcript_ids'] == 'G101.1']
        self.assertTrue(all(row['coverage_pass'] == 'False' for row in bad))
        exact = next(row for row in models if row['transcript_ids'] == 'G100.1')
        self.assertEqual(exact['n_required_bases'], '75')
        summary = {row['metric']: row['value'] for row in self.table('.summary.tsv')}
        self.assertEqual(summary['n_regions_excluded_coordinate_convention'], '1')
        self.assertEqual(summary['n_final_regions'], '5')
        with (self.directory/'reference.dapars2.bed').open() as source:
            native = list(csv.reader(source, delimiter='\t'))
        with (self.directory/'reference.full_exon.igv.bed').open() as source:
            igv = list(csv.reader(source, delimiter='\t'))
        for a, b in zip(native, igv):
            self.assertEqual(int(a[1]), int(b[1])+1)
            self.assertEqual(a[2:], b[2:])

    def test_parallel_and_serial_outputs_agree(self):
        self.run_cli(prefix='serial', workers=1)
        self.run_cli(prefix='parallel', workers=2)
        for suffix in ('.dapars2.bed', '.full_exon.igv.bed', '.regions.tsv.gz', '.terminal_models.tsv.gz', '.summary.tsv'):
            a, b = self.directory/('serial'+suffix), self.directory/('parallel'+suffix)
            opener = gzip.open if a.suffix == '.gz' else open
            with opener(a, 'rb') as left, opener(b, 'rb') as right:
                self.assertEqual(left.read(), right.read())

    def test_mapped_depth_normalization_and_ceiling_support(self):
        # The 100M-alignment sample needs raw depth >=4, so depth 3 must fail.
        self.write_bigwig(self.directory/'s1.bw', 3)
        self.run_cli()
        self.assertEqual(len(self.table('.regions.tsv.gz')), 0)
        self.assertTrue(all(row['n_pass'] in {'0', '1'} for row in self.table('.terminal_models.tsv.gz')))

    def test_prefix_integration_matches_per_base_counting(self):
        runs = np.array([[10, 15], [20, 30], [40, 42]], dtype=np.int64)
        coverage = np.zeros(60, dtype=bool)
        for start, end in runs:
            coverage[start:end] = True
        points = np.arange(61)
        np.testing.assert_array_equal(refprep.covered_to(points, runs), np.r_[0, coverage.cumsum()])
        self.assertEqual(refprep.required_count(100, 0.07), 7)
        self.assertEqual(refprep.required_count(6, 0.2), 2)

    def test_gene_biotype_without_gene_records(self):
        with gzip.open(self.gtf, 'rt') as source:
            content = source.read().replace('gene_type ', 'gene_biotype ')
        with gzip.open(self.gtf, 'wt') as output:
            output.write(content)
        self.run_cli(extra=('--chromosomes', 'chr1'))
        self.assertEqual(len(self.table('.regions.tsv.gz')), 5)
        self.assertTrue(all(row['chrom'] == 'chr1' for row in self.table('.terminal_models.tsv.gz')))

    def test_gene_records_can_supply_types_and_symbols_after_exons(self):
        with gzip.open(self.gtf, 'rt') as source:
            content = source.read()
        content = re.sub(r'gene_type "[^"]*";\s*', '', content)
        content = re.sub(r'gene_name "[^"]*";\s*', '', content)
        for number in range(1, 8):
            chrom = 'chr2' if number == 6 else 'chr1'
            strand = '-' if number in (2, 6) else '+'
            gene_type = 'lncRNA' if number == 4 else 'protein_coding'
            content += (f'{chrom}\ttest\tgene\t1\t1000\t.\t{strand}\t.\t'
                        f'gene_id "ENSG{number:011d}.2_7"; gene_name "GENE{number}"; gene_type "{gene_type}";\n')
        with gzip.open(self.gtf, 'wt') as output:
            output.write(content)
        self.run_cli()
        regions = self.table('.regions.tsv.gz')
        self.assertEqual(len(regions), 5)
        self.assertTrue(all(row['gene_symbol'].startswith('GENE') for row in regions))

    def test_reserved_transcript_separator_is_not_silently_reinterpreted(self):
        with gzip.open(self.gtf, 'rt') as source:
            content = source.read().replace('A1.1', 'A1,fragment.1')
        with gzip.open(self.gtf, 'wt') as output:
            output.write(content)
        result = self.run_cli(success=False)
        self.assertIn('transcript IDs cannot contain commas', result.stderr)

    def test_duplicate_samples_do_not_inflate_support(self):
        with self.samples.open('a') as output:
            output.write('duplicate\ts0.bw\t50000000\n')
        result = self.run_cli(success=False)
        self.assertIn('must be unique', result.stderr)

    def test_missing_chromosome_and_negative_signal_are_errors(self):
        self.write_bigwig(self.directory/'s0.bw', 2, chromosomes=('chr1',))
        result = self.run_cli(success=False)
        self.assertIn('chromosome chr2', result.stderr)
        self.write_bigwig(self.directory/'s0.bw', -2)
        result = self.run_cli(success=False)
        self.assertIn('nonnegative raw read coverage', result.stderr)

    def test_gene_id_matching_does_not_strip_arbitrary_dots(self):
        self.assertEqual(refprep.stable_gene_id('ENSG00000000001.2_7'), 'ENSG00000000001')
        self.assertEqual(refprep.stable_gene_id('ENSMUSG00000000001.3'), 'ENSMUSG00000000001')
        self.assertEqual(refprep.stable_gene_id('custom.gene.2'), 'custom.gene.2')


if __name__ == '__main__':
    unittest.main()
