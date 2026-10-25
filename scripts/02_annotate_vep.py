"""Annotate candidates with Ensembl VEP REST (CADD, SpliceAI, AlphaMissense, REVEL, gnomAD, ClinVar)
and apply the recessive rare-variant filter. Only coordinates are sent (no identifiers)."""
import sys, csv, json, urllib.request
rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
rows = [r for r in rows if 'PASS' in r['filter']]
body = json.dumps({'variants': [f"{r['chrom']} {r['pos']} . {r['ref']} {r['alt']} . . ." for r in rows]}).encode()
url = 'https://rest.ensembl.org/vep/homo_sapiens/region?CADD=1&SpliceAI=1&AlphaMissense=1&REVEL=1&canonical=1&hgvs=1&pick=1'
req = urllib.request.Request(url, body, {'Content-Type': 'application/json', 'Accept': 'application/json'})
ann = json.load(urllib.request.urlopen(req))
LOF = ('stop_gained', 'frameshift', 'splice_donor', 'splice_acceptor', 'start_lost', 'missense', 'inframe')
print('gene\tvariant\tconsequence\thgvsp\tcadd\talphamissense\tspliceai_max\tgnomad_af\tclinvar\tpass_filter')
for x in ann:
    t = (x.get('transcript_consequences') or [{}])[0]
    cv = next((c for c in x.get('colocated_variants', []) if c.get('frequencies')), {})
    af = list(cv['frequencies'].values())[0].get('gnomadg') if cv else None
    sai = t.get('spliceai') or {}
    smax = max([v for k, v in sai.items() if k.startswith('DS')], default=0)
    rare = af is None or af < 1e-3
    keep = rare and (any(k in x['most_severe_consequence'] for k in LOF) or smax >= 0.2)
    print(t.get('gene_symbol'), x['input'].split()[0] + ':' + x['input'].split()[1], x['most_severe_consequence'],
          t.get('hgvsp'), t.get('cadd_phred'), (t.get('alphamissense') or {}).get('am_pathogenicity'), smax, af,
          '/'.join(cv.get('clin_sig', [])), keep, sep='\t')
