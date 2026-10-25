"""Extract non-reference calls in MVA genes (GRCh38) from the proband VCF.
Usage: python 01_extract_candidates.py proband.vcf.gz > results/candidates.tsv
Patient data is NOT included in this repository (controlled access)."""
import sys, pysam
GENES = {'BUB1B': ('15', 40150000, 40230000), 'TRIP13': ('5', 880000, 925000), 'CEP57': ('11', 95780000, 95840000)}
v = pysam.VariantFile(sys.argv[1])
print('gene\tchrom\tpos\tref\talt\tgt\tfilter\tad\tdp')
for g, (c, s, e) in GENES.items():
    for r in v.fetch(c, s, e):
        gt = r.samples[0]['GT']
        if gt in ((0, 0), (None, None)):
            continue
        print(g, r.chrom, r.pos, r.ref, ','.join(r.alts), '/'.join(map(str, gt)), ';'.join(r.filter.keys()),
              ','.join(map(str, r.samples[0].get('AD') or [])), r.samples[0].get('DP'), sep='\t')
