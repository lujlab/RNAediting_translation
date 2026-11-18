#!/usr/bin/env python3
"""Estimate half-lives from published expression and annotate editing in FlyBase r6.12.

The embryonic editing table uses dm6 coordinates. Select the longest CDS from
transcripts containing a detected site for edited genes and from all CDS
transcripts for controls. Break ties by mature transcript length, then ID.
Pool positive editing intervals for sites in the selected CDS.
"""
import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, required=True)
parser.add_argument('--gtf', type=Path, required=True)
parser.add_argument('--genome-fasta', type=Path, required=True)
parser.add_argument('--transcript-fasta', type=Path, required=True)
parser.add_argument('--editing-file', type=Path, required=True)
parser.add_argument('--half-life-file', type=Path, required=True)
args = parser.parse_args()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
GTF = args.gtf.resolve(strict=True)
GENOME = args.genome_fasta.resolve(strict=True)
TX_FASTA = args.transcript_fasta.resolve(strict=True)
LONG = args.editing_file.resolve(strict=True)
HALF_LIFE = args.half_life_file.resolve(strict=True)


def read(path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0])
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t',
                           lineterminator='\n', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def merge(blocks):
    result = []
    for lo, hi in sorted(set(blocks)):
        if result and lo <= result[-1][1] + 1:
            result[-1] = (result[-1][0], max(hi, result[-1][1]))
        else:
            result.append((lo, hi))
    return result


def fasta(path):
    name, chunks = None, []
    with path.open(encoding='utf-8') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks).upper()


# The CG9243 microarray measurement is assigned once to the edited ORF CG43346.
alias_rows = []
alias_genes = defaultdict(set)
legacy_cg_to_gene = {}
with TX_FASTA.open(encoding='utf-8') as f:
    for header in f:
        if not header.startswith('>'):
            continue
        tid = header[1:].split()[0]
        parent = re.search(r'parent=(FBgn\d+)', header)
        name = re.search(r'name=([^;]+)', header)
        if not parent:
            continue
        header_aliases = re.findall(r'FlyBase_Annotation_IDs:(CG\d+)', header)
        if header_aliases:
            legacy_cg_to_gene.setdefault(header_aliases[-1], parent[1])
        for alias in sorted(set(header_aliases)):
            alias_genes[alias].add(parent[1])
            alias_rows.append(dict(annotation_id=alias, annotation_FBgn=parent[1],
                transcript_id=tid, transcript_symbol=name[1] if name else ''))
write(OUT/'r6_12_all_historical_CG_aliases.tsv', alias_rows)
assert alias_genes['CG9243'] == {'FBgn0263050', 'FBgn0263051'}
assert any(r['annotation_id']=='CG9243' and r['annotation_FBgn']=='FBgn0263051'
           and r['transcript_id']=='FBtr0306975' for r in alias_rows)

# Use the last CG identifier in each transcript header for the source-gene mapping.
gene_symbols = {}
legacy_exons = {}
with GTF.open(encoding='utf-8') as f:
    for line in f:
        a = line.rstrip().split('\t')
        if len(a) != 9 or a[2] not in ('gene', 'exon'):
            continue
        attrs = dict(re.findall(r'(\w+)\s+"([^"]*)"', a[8]))
        gid, tid = attrs.get('gene_id'), attrs.get('transcript_id')
        if gid and attrs.get('gene_symbol'):
            gene_symbols.setdefault(gid, attrs['gene_symbol'])
        if gid and tid and a[2] == 'exon':
            item = legacy_exons.setdefault(tid, [gid, 0])
            assert item[0] == gid
            item[1] += int(a[4]) - int(a[3]) + 1
symbol_to_gene = {}
for gid, symbol in sorted(gene_symbols.items()):
    symbol_to_gene.setdefault(symbol, gid)
legacy_transcripts = {}
for tid, (gid, length) in legacy_exons.items():
    candidate = (length, tid)
    old = legacy_transcripts.get(gid)
    if old is None or (-length, tid) < (-old[0], old[1]):
        legacy_transcripts[gid] = candidate

def clean_id(value):
    return re.sub(r'\.\d+$', '', value.strip())

def numeric(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float('nan')

def finite_median(values):
    values = [v for v in values if math.isfinite(v)]
    return statistics.median(values) if values else ''

raw_half_life = read(HALF_LIFE)
assert raw_half_life, 'Published half-life table is empty'
series_cols = {label: [k for k in raw_half_life[0] if re.match('^'+label+r'[._]', k)]
               for label in ('U2', 'U3', 'E2', 'E3')}
assert all(series_cols.values()), 'Expected normalized U2/U3/E2/E3 replicate columns'
probe_rows = []
probes_by_gene = defaultdict(list)
for r in raw_half_life:
    decay_class = r['class.I.to.V'].strip().upper()
    if decay_class not in ('I', 'II', 'III', 'IV', 'V'):
        continue
    source_cg = clean_id(r['CG.gene.ID'])
    gid = source_cg if source_cg.startswith('FBgn') else legacy_cg_to_gene.get(source_cg)
    if not gid:
        gid = symbol_to_gene.get(source_cg) or symbol_to_gene.get(clean_id(r.get('Gene.symbol','')))
    if not gid:
        continue
    means = {}
    for label, columns in series_cols.items():
        values = [numeric(r[k]) for k in columns]
        values = [v for v in values if math.isfinite(v)]
        means[label] = sum(values)/len(values) if values else float('nan')
    if decay_class in ('II', 'III'):
        difference = means['U2'] - means['U3']
    elif decay_class in ('IV', 'V'):
        difference = means['E2'] - means['E3']
    else:
        difference = float('nan')
    half = 60.0/difference if math.isfinite(difference) and difference > 0 else float('nan')
    baseline = means['U2'] if decay_class in ('II', 'III') else means['E2']
    expr = 2.0**baseline if math.isfinite(baseline) else float('nan')
    p = dict(source_gene=source_cg, FBgn=gid, decay_class=decay_class,
             probe_set=r.get('Affy.probe.set.ID',''), half_life_min=half,
             published_expression=expr, log2_baseline_expression=baseline)
    probe_rows.append(p)
    probes_by_gene[gid].append(p)
conflicts = [dict(FBgn=g, n_probe_sets=len(rows),
    n_decay_classes=len({r['decay_class'] for r in rows}),
    classes=';'.join(sorted({r['decay_class'] for r in rows})))
    for g, rows in sorted(probes_by_gene.items()) if len({r['decay_class'] for r in rows}) > 1]
conflict_ids = {r['FBgn'] for r in conflicts}
write(OUT/'2010_conflicting_probe_classifications.tsv', conflicts,
      ['FBgn','n_probe_sets','n_decay_classes','classes'])
write(OUT/'2010_probe_level_reestimated.tsv', probe_rows)
raw_editing_rows = read(LONG)
original_long = [r for r in raw_editing_rows if float(r['ratio']) > 0]
assert original_long, 'No positive-editing stage records'
edited_source_ids = {clean_id(r['gene']) for r in original_long}
source_meta = []
for gid, probes in sorted(probes_by_gene.items()):
    if gid in conflict_ids:
        continue
    assert gid in legacy_transcripts, ('No source representative transcript', gid)
    decay_class = probes[0]['decay_class']
    length, tid = legacy_transcripts[gid]
    source_meta.append(dict(FBgn=gid, gene=probes[0]['source_gene'],
        gene_symbol=gene_symbols.get(gid,''), decay_class=decay_class,
        half_life_min=finite_median([p['half_life_min'] for p in probes]),
        published_expression=finite_median([p['published_expression'] for p in probes]),
        log2_baseline_expression=finite_median([p['log2_baseline_expression'] for p in probes]),
        expression_source='E-MEXP-2580 processed U2 mean (class-specific decay-series baseline)'
            if decay_class in ('II','III') else 'E-MEXP-2580 processed E2 mean (embryo-series baseline)',
        n_probe_sets=len(probes), source_probe_sets=';'.join(p['probe_set'] for p in probes),
        transcript_id=tid, transcript_length=length,
        editing_binary_status='Edited' if gid in edited_source_ids else 'No_detected_editing'))
write(OUT/'2010_source_gene_measurements.tsv', source_meta)

gene_overrides = {'FBgn0263050': 'FBgn0263051'}
assert not any(r['FBgn']=='FBgn0263051' for r in source_meta), 'Would duplicate a measurement'
meta = []
for source in source_meta:
    row = dict(source)
    target = gene_overrides.get(source['FBgn'], source['FBgn'])
    row.update(source_FBgn=source['FBgn'], previous_assigned_FBgn=source['FBgn'],
        measurement_gene=source['gene'], source_gene_symbol=source['gene_symbol'],
        annotation_FBgn=target, FBgn=target,
        measurement_scope='historical_gene_level',
        annotation_mapping_note='CG9243 edited ORF resolved to CG43346; historical measurement retained once'
            if source['FBgn'] in gene_overrides else 'unchanged')
    if source['FBgn'] in gene_overrides:
        assert source['gene']=='CG9243'
        row['gene_symbol']='CG43346'
        row['measurement_probe_set']='1632466_at'
    else:
        row['measurement_probe_set']=''
    meta.append(row)
meta_by_gene = {r['FBgn']: r for r in meta}
assert len(meta_by_gene) == len(meta)
long = []
for source in original_long:
    row = dict(source)
    row['source_FBgn'] = source['gene']
    row['gene'] = gene_overrides.get(source['gene'], source['gene'])
    if source['gene'] in gene_overrides:
        assert source['gene_symbol_2011']=='CG9243'
        row['gene_symbol_r6_12']='CG43346'
    long.append(row)
write(OUT/'resolved_embryo_CDS_editing_sites_long.tsv', long)
mapping_audit = []
for r in meta:
    if r['editing_binary_status']!='Edited':
        continue
    names = {r['measurement_gene'], r['source_gene_symbol']}
    names.update(s['gene_symbol_2011'] for s in long if s['gene']==r['FBgn'])
    all_ids = set().union(*(alias_genes.get(name,set()) for name in names))
    mapping_audit.append(dict(measurement_gene=r['measurement_gene'], source_FBgn=r['source_FBgn'],
        annotation_FBgn=r['annotation_FBgn'], historical_names=';'.join(sorted(names)),
        historical_alias_FBgns=';'.join(sorted(all_ids)), n_alias_FBgns=len(all_ids),
        mapping_note=r['annotation_mapping_note']))
write(OUT/'edited_gene_historical_alias_audit.tsv', mapping_audit)
dedup = {}
for r in long:
    key = (r['gene'], r['dmel_chr'], r['dmel_pos'], r['lib'])
    if key not in dedup or int(r['depth']) > int(dedup[key]['depth']):
        dedup[key] = r
pooled = {}
for r in dedup.values():
    key = (r['gene'], r['dmel_chr'], r['dmel_pos'])
    p = pooled.setdefault(key, dict(gene=r['gene'], dmel_chr=r['dmel_chr'],
                                   dmel_pos=int(r['dmel_pos']), altcount=0, depth=0))
    p['altcount'] += int(r['altcount'])
    p['depth'] += int(r['depth'])

print('Reading matching r6.12 genomic DNA', flush=True)
genome = dict(fasta(GENOME))
sites = defaultdict(list)
for p in pooled.values():
    if p['gene'] not in meta_by_gene:
        continue
    assert 0 <= p['altcount'] <= p['depth'] and p['depth'] > 0
    p['site_editing_fraction'] = p['altcount'] / p['depth']
    chrom = p['dmel_chr'][3:] if p['dmel_chr'].startswith('chr') else p['dmel_chr']
    base = genome[chrom][p['dmel_pos'] - 1]
    assert base in 'AT', (p['gene'], chrom, p['dmel_pos'], base)
    p.update(chromosome=chrom, strand='+' if base == 'A' else '-')
    sites[p['gene']].append(p)

models = {}
with GTF.open(encoding='utf-8') as f:
    for line in f:
        if line.startswith('#'):
            continue
        a = line.rstrip().split('\t')
        if len(a) != 9 or a[2] not in ('exon', 'CDS', '5UTR', '3UTR'):
            continue
        attrs = dict(re.findall(r'(\w+)\s+"([^"]*)"', a[8]))
        gene, tid = attrs.get('gene_id'), attrs.get('transcript_id')
        if gene not in meta_by_gene or not tid:
            continue
        m = models.setdefault(tid, dict(gene=gene, transcript_id=tid,
            transcript_symbol=attrs.get('transcript_symbol', ''), chrom=a[0],
            strand=a[6], regions=defaultdict(list)))
        assert (m['gene'], m['chrom'], m['strand']) == (gene, a[0], a[6])
        m['regions'][a[2]].append((int(a[3]), int(a[4])))

bygene = defaultdict(list)
candidates = []
for m in models.values():
    for region, blocks in m['regions'].items():
        m['regions'][region] = merge(blocks)
    m['transcript_length'] = sum(hi-lo+1 for lo, hi in m['regions'].get('exon', []))
    m['CDS_length'] = sum(hi-lo+1 for lo, hi in m['regions'].get('CDS', []))
    if not m['transcript_length'] or not m['CDS_length']:
        continue
    for lo, hi in m['regions']['CDS']:
        assert any(a <= lo <= hi <= b for a, b in m['regions']['exon']), m['transcript_id']
    m['covered_sites'] = [s for s in sites[m['gene']]
        if s['chromosome'] == m['chrom'] and s['strand'] == m['strand']
        and any(lo <= s['dmel_pos'] <= hi for lo, hi in m['regions']['CDS'])]
    m['n_sites'] = len(m['covered_sites'])
    bygene[m['gene']].append(m)
    candidates.append(dict(FBgn=m['gene'], transcript_id=m['transcript_id'],
        transcript_symbol=m['transcript_symbol'], chromosome=m['chrom'], strand=m['strand'],
        transcript_length=m['transcript_length'], CDS_length=m['CDS_length'],
        n_CDS_sites=m['n_sites']))

complement = str.maketrans('ACGTN', 'TGCAN')
def sequence(m, region):
    s = ''.join(genome[m['chrom']][lo-1:hi] for lo, hi in m['regions'].get(region, []))
    return s.translate(complement)[::-1] if m['strand'] == '-' else s


features, excluded, site_audit, region_rows = [], [], [], []
selected, exon_sequences = {}, {}
with (OUT / 'selected_CDS_sequences.fasta').open('w', encoding='utf-8') as ff:
    for gene, old in meta_by_gene.items():
        was_edited = old['editing_binary_status'] == 'Edited'
        eligible = bygene.get(gene, [])
        if was_edited:
            eligible = [m for m in eligible if m['n_sites'] > 0]
        m = min(eligible, key=lambda x: (-x['CDS_length'], -x['transcript_length'], x['transcript_id'])) if eligible else None
        for s in sites[gene]:
            site_audit.append(dict(FBgn=gene, source_FBgn=old['source_FBgn'],
                measurement_gene=old['measurement_gene'], annotation_FBgn=old['annotation_FBgn'],
                gene_symbol=old['gene_symbol'],
                transcript_id=m['transcript_id'] if m else '',
                dmel_chr=s['chromosome'], dmel_pos=s['dmel_pos'], strand=s['strand'],
                altcount=s['altcount'], depth=s['depth'], editing_fraction=s['site_editing_fraction'],
                within_any_gene_CDS=any(s in x['covered_sites'] for x in bygene.get(gene, [])),
                included_in_selected_CDS=bool(m and s in m['covered_sites'])))
        if m is None:
            excluded.append(dict(FBgn=gene, gene=old['gene'], gene_symbol=old['gene_symbol'],
                previous_status=old['editing_binary_status'],
                reason='No r6.12 CDS transcript covers a detected input site' if was_edited
                       else 'No r6.12 CDS-bearing transcript'))
            continue
        selected[m['transcript_id']] = m
        exon_seq, cds_seq = sequence(m, 'exon'), sequence(m, 'CDS')
        exon_sequences[m['transcript_id']] = exon_seq
        assert len(exon_seq) == m['transcript_length'] and len(cds_seq) == m['CDS_length']
        assert cds_seq.count('A') > 0
        levels = [s['site_editing_fraction'] for s in m['covered_sites']] if was_edited else []
        alt = sum(s['altcount'] for s in m['covered_sites']) if was_edited else 0
        depth = sum(s['depth'] for s in m['covered_sites']) if was_edited else 0
        row = dict(old)
        row.update(transcript_id=m['transcript_id'], transcript_symbol=m['transcript_symbol'],
            annotation_source='FlyBase_r6.12_GTF', reference_assembly='dm6',
            transcript_length=m['transcript_length'], CDS_length=m['CDS_length'],
            UTR5_length=sum(hi-lo+1 for lo, hi in m['regions'].get('5UTR', [])),
            UTR3_length=sum(hi-lo+1 for lo, hi in m['regions'].get('3UTR', [])),
            exon_number=len(m['regions']['exon']), intron_number=len(m['regions']['exon'])-1,
            total_intron_length=m['regions']['exon'][-1][1]-m['regions']['exon'][0][0]+1-m['transcript_length'],
            GC_content=(exon_seq.count('G')+exon_seq.count('C'))/len(exon_seq),
            A_content=exon_seq.count('A')/len(exon_seq), editable_A_count=exon_seq.count('A'),
            CDS_A_count=cds_seq.count('A'), editing_site_count=len(levels), sum_editing=sum(levels),
            edited_allele_count=alt, total_allele_count=depth,
            mean_editing=statistics.mean(levels) if levels else '',
            median_editing=statistics.median(levels) if levels else '',
            max_editing=max(levels) if levels else '',
            mean_depth=statistics.mean(s['depth'] for s in m['covered_sites']) if levels else '',
            gene_editing_burden=alt/depth if depth else '',
            original_transcript_id=old['transcript_id'], original_transcript_length=old['transcript_length'],
            n_input_gene_CDS_sites=len(sites[gene]), all_input_sites_covered=len(levels)==len(sites[gene]))
        features.append(row)
        ff.write('>{} FBgn={} source=FlyBase_r6.12_GTF\n{}\n'.format(m['transcript_id'], gene, cds_seq))
        for region, blocks in m['regions'].items():
            for lo, hi in blocks:
                region_rows.append(dict(FBgn=gene, transcript_id=m['transcript_id'],
                    chromosome=m['chrom'], strand=m['strand'], region=region, start=lo, end=hi))

tx_checks = []
for tid, seq in fasta(TX_FASTA):
    if tid in exon_sequences:
        tx_checks.append(dict(transcript_id=tid, genomic_length=len(exon_sequences[tid]),
                              fasta_length=len(seq), exact_sequence_match=seq==exon_sequences[tid]))
assert len(tx_checks) == len(selected)
write(OUT/'transcript_FASTA_QC.tsv', tx_checks)
assert all(r['exact_sequence_match'] for r in tx_checks), 'r6.12 transcript/genome sequence mismatch'
for r in candidates:
    r['selected'] = r['transcript_id'] in selected
write(OUT/'selected_gene_features.tsv', features)
write(OUT/'transcript_candidates.tsv', candidates)
write(OUT/'selected_transcript_site_audit.tsv', site_audit)
write(OUT/'selected_transcript_regions.tsv', region_rows)
write(OUT/'excluded_genes.tsv', excluded, ['FBgn','gene','gene_symbol','previous_status','reason'])
with GTF.open(encoding='utf-8') as src, (OUT/'selected_transcripts_r6.12.gtf').open('w', encoding='utf-8') as dst:
    dst.write('# Native FlyBase r6.12 coordinates; selected representative transcripts.\n')
    for line in src:
        match = re.search(r'transcript_id "([^"]+)"', line)
        if match and match[1] in selected:
            dst.write(line)

def valid_half_life(r):
    return r['decay_class'] in ('II','III','IV','V') and r['half_life_min'] != '' and math.isfinite(float(r['half_life_min'])) and float(r['half_life_min']) > 0
orig_e = [r for r in meta if r['editing_binary_status']=='Edited' and valid_half_life(r)]
new_e = [r for r in features if r['editing_binary_status']=='Edited' and valid_half_life(r)]
e_genes = {r['FBgn'] for r in orig_e}
write(OUT/'panel_E_excluded_input_sites.tsv', [r for r in site_audit if r['FBgn'] in e_genes and not r['included_in_selected_CDS']])
summary = dict(annotation='FlyBase r6.12 GTF; dm6 genomic DNA',
    mapped_class_gene_count=len(probes_by_gene), conflicting_class_gene_count=len(conflicts),
    measurement_source='Published normalized microarray signals; half-lives re-estimated from replicate means',
    mapping_correction='Historical CG9243 retained once; edited CDS annotation CG43346 / FBgn0263051',
    original_gene_count=len(meta), retained_gene_count=len(features),
    original_edited=sum(r['editing_binary_status']=='Edited' for r in meta),
    retained_edited=sum(r['editing_binary_status']=='Edited' for r in features),
    original_E_genes=len(orig_e), retained_E_genes=len(new_e),
    original_E_sites=sum(len(sites[r['FBgn']]) for r in orig_e),
    retained_E_sites=sum(r['editing_site_count'] for r in new_e),
    excluded_genes=len(excluded), excluded_edited=[r['gene_symbol'] for r in excluded],
    long_table_rows=len(long), input_editing_rows=len(raw_editing_rows),
    excluded_nonpositive_editing_rows=len(raw_editing_rows)-len(original_long),
    deduplicated_rows=len(dedup),
    input_library_count=len({r['lib'] for r in long}),
    zero_editing_rows=sum(float(r['ratio'])<=0 for r in long),
    matched_transcript_FASTA_sequences=len(tx_checks),
    selected_CDS_N_bases=sum(s.count('N') for tid,s in fasta(OUT/'selected_CDS_sequences.fasta')))
(OUT/'annotation_QC.json').write_text(json.dumps(summary,indent=2))
manifest = []
for p in [GTF,GENOME,TX_FASTA,HALF_LIFE,LONG]:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            h.update(chunk)
    manifest.append(dict(path=str(p),bytes=p.stat().st_size,sha256=h.hexdigest()))
(OUT/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(summary,indent=2),flush=True)
