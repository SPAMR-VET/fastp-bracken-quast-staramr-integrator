###################################################
#
# Includes Fastp, Bracken, and Quast results into STARAMR spreadsheet output
# Input files are tooldistillator outputs from the Abromics Workflow, in JSON or EXCEL format.
# Perform quality check at the same time.
#
# Created: 2025-01-29
# Last Edited: 2025-05-14
#
# Contributors:
#
# Elena Lavinia Diaconu (elena.diaconu@izslt.it)
# Maciej Kochanowski (maciej.kochanowski@piwet.pulawy.pl)
# Sebastien Leclercq (sebastien.leclercq@inrae.fr)
#
# developed as part of the SPAMR-VET Consortium (EUPAHW)
###################################################

import os
import re
import json
import argparse
import numpy as np
import pandas as pd
from openpyxl import load_workbook

def extract_sample_id(file_path):
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        if isinstance(data, list):
            fastp_section = next((i for i in data if i.get('analysis_software_name') == 'fastp'), None)
            if fastp_section:
                fastp_report = next((r for r in fastp_section.get('results', []) if r.get('name') == 'fastp_report'), None)
                if fastp_report and fastp_report.get('content'):
                    command_str = fastp_report['content'][0].get('command', '')
                    match = re.search(r'fastp report for ([^\s]+)\.fastqsanger\.gz', command_str)
                    if match:
                        return match.group(1)
        
        if isinstance(data, dict) and data.get('analysis_software_name') == 'quast':
            for r in data.get('results', []):
                if r.get('name') == 'quast_report' and isinstance(r.get('content'), list) and r['content']:
                    assembly = r['content'][0].get('Assembly', None)
                    if assembly:
                        return assembly
        
        if isinstance(data, list):
            quast_section = next((i for i in data if i.get('analysis_software_name') == 'quast'), None)
            if quast_section and 'results' in quast_section:
                for r in quast_section['results']:
                    if r.get('name') == 'quast_report' and isinstance(r.get('content'), list) and r['content']:
                        assembly = r['content'][0].get('Assembly', None)
                        if assembly:
                            return assembly
    
    except Exception as e:
        pass
    
    base = os.path.basename(file_path)
    name = re.sub(r'\.(json|csv|tsv|fasta|fa|xlsx|txt|dat)$', '', base, flags=re.IGNORECASE)
    for suffix in ['_FastpKraken', '_quast', '_fastp', '_kraken', '_bracken', '_recentrifuge', '_report', '-1']:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break
    return name.strip('_- ')

def safe_float_convert(v, default=np.nan):
    try:
        if isinstance(v, str):
            v = v.replace(',', '')
        return float(v)
    except:
        return default

def safe_int_convert(v, default=np.nan):
    try:
        if isinstance(v, str):
            v = v.replace(',', '')
        return int(float(v))
    except:
        return default

def parse_fastp_json(data):
    sec = next((i for i in data if i.get('analysis_software_name') == 'fastp'), None)
    if not sec:
        raise ValueError("fastp section missing")
    rep = next((r for r in sec.get('results', []) if r.get('name') == 'fastp_report'), None)
    if not rep or not rep.get('content'):
        raise ValueError("fastp_report missing")
    c = rep['content'][0]
    pre = c.get('summary', {}).get('before_filtering', {})
    post = c.get('summary', {}).get('after_filtering', {})
    filt = c.get('filtering_result', {})
    dup = c.get('duplication', {})
    return {
        'total_reads': pre.get('total_reads'),
        'total_bases': pre.get('total_bases'),
        'q20_rate': pre.get('q20_rate'),
        'q30_rate': pre.get('q30_rate'),
        'gc_content': pre.get('gc_content'),
        'duplication_rate': dup.get('rate'),
        'passed_filter_reads': filt.get('passed_filter_reads'),
        'low_quality_reads': filt.get('low_quality_reads'),
        'total_reads_after': post.get('total_reads'),
        'ratio_reads_after': post.get('total_reads')/pre.get('total_reads') if pre.get('total_reads') else np.nan,
        'total_bases_after': post.get('total_bases'),
        'q30_rate_after': post.get('q30_rate'),
        'gc_content_after': post.get('gc_content')
    }

def extract_top_taxa_from_bracken(data):
    sec = next((i for i in data if i.get('analysis_software_name') == 'bracken'), None)
    result = {f'{p}_{i}': np.nan for i in range(1,6) for p in ['taxon','count','proportion']}
    if not sec:
        return result
    rep = next((r for r in sec.get('results', []) if r.get('name') == 'bracken_report'), None)
    content = rep.get('content', []) if rep else []
    taxa = []
    for item in content:
        taxa.append({
            'name': item['name'],
            'count': item.get('new_est_reads'),
            'score': item.get('fraction_total_reads')
        })
    taxa = sorted(taxa, key=lambda x: x['count'] or 0, reverse=True)[:5]
    for idx, t in enumerate(taxa, 1):
        result[f'taxon_{idx}'] = t['name']
        result[f'count_{idx}'] = safe_int_convert(t['count'])
        result[f'proportion_{idx}'] = safe_float_convert(t['score'])*100
    return result

def parse_quast_json(path):
    try:
        with open(path, 'r') as f:
            data = json.load(f)
        sec = None
        if isinstance(data, list):
            sec = next((i for i in data if i.get('analysis_software_name') == 'quast'), None)
        elif isinstance(data, dict) and data.get('analysis_software_name') == 'quast':
            sec = data
        if sec and 'results' in sec:
            for r in sec['results']:
                if r.get('name') == 'quast_report' and isinstance(r.get('content'), list) and r['content']:
                    raw = r['content'][0]
                    break
            else:
                raise ValueError("quast_report content missing")
        else:
            raise ValueError("quast section missing")
        cleaned = {}
        for k, v in raw.items():
            key = str(k).strip()
            key = key.replace(' ', '_').replace('#', 'num').replace("(%", 'percent').replace(")", '').replace("(", '')
            key = re.sub(r'contigs_>=_?(\d+)_?bp', r'num_contigs_\1bp', key)
            if key == 'contigs':
                key = 'num_contigs_500bp'
            if key == 'GC':
                key = 'GC_percent'
            cleaned[key] = v
        for e in ['N50','Total_length','Largest_contig','num_contigs_500bp','GC_percent']:
            if e not in cleaned:
                cleaned[e] = np.nan
        return cleaned
    except Exception:
        try:
            df = pd.read_csv(path, sep='\t', engine='python')
        except:
            df = pd.read_csv(path, sep=',', engine='python')
        if df.shape[1] == 2 and df.shape[0] > 1:
            d = pd.Series(df.iloc[:,1].values, index=df.iloc[:,0]).to_dict()
        else:
            d = df.iloc[0].to_dict()
        cleaned = {}
        for k,v in d.items():
            key = str(k).strip().replace(' ', '_').replace('#','num').replace('(%)','percent')
            key = re.sub(r'contigs_>=_?(\d+)_?bp', r'num_contigs_\1bp', key)
            if key == 'GC':
                key = 'GC_percent'
            cleaned[key] = v
        for e in ['N50','Total_length','Largest_contig','num_contigs_500bp','GC_percent']:
            if e not in cleaned:
                cleaned[e] = np.nan
        return cleaned

def calculate_quality_metrics(fastp, bracken, quast):
    thresholds = {
        'total reads after filtering': {'warning': 500000, 'fail': 400000},
        'ratio of filtered reads': {'warning': 0.8, 'fail': 0.6},
        'Bracken main taxon proportion': {'warning': 90, 'fail': 80},
        'reads gc content': {'warning': (0.3,0.7), 'fail': (0.25,0.75)},
        'N50': {'warning': 35000, 'fail': 30000},
        'largest contig': {'warning': 100000, 'fail': 50000},
        'number of contigs >500 bp': {'warning': 400, 'fail': 500},
        'number of contigs >0 bp': {'warning': 500, 'fail': 700},
    }
    fb = []
    result = 'Passed'
    def chk(val, name, mode='low'):
        nonlocal result
        if pd.isna(val):
            fb.append(f"Missing {name}")
            return
        w = thresholds[name]['warning']
        f = thresholds[name]['fail']
        if mode == 'low':
            if val < f:
                result, fb_msg = 'Fail', f"Low {name} (<{f})"
            elif val < w:
                if result != 'Fail':
                    result = 'Warning'
                fb_msg = f"Low {name} (<{w})"
            else:
                return
        elif mode == 'high':
            if val > f:
                result, fb_msg = 'Fail', f"High {name} (>{f})"
            elif val > w:
                if result != 'Fail':
                    result = 'Warning'
                fb_msg = f"High {name} (>{w})"
            else:
                return
        elif mode == 'range':
            if not (f[0] <= val <= f[1]):
                result, fb_msg = 'Fail', f"Abnormal {name}"
            elif not (w[0] <= val <= w[1]):
                if result != 'Fail':
                    result = 'Warning'
                fb_msg = f"Borderline {name}"
            else:
                return
        fb.append(fb_msg)
    
    chk(safe_int_convert(fastp.get('total_reads_after')), 'total reads after filtering', 'low')
    chk(safe_float_convert(fastp.get('ratio_reads_after')), 'ratio of filtered reads', 'low')
    chk(safe_float_convert(bracken.get('proportion_1')), 'Bracken main taxon proportion', 'low')
    chk(safe_float_convert(fastp.get('gc_content')), 'reads gc content', 'range')
    chk(safe_int_convert(quast.get('N50')), 'N50', 'low')
    chk(safe_int_convert(quast.get('Largest_contig')), 'largest contig','low')
    chk(safe_int_convert(quast.get('num_contigs_500bp')), 'number of contigs >500 bp', 'high')
    chk(safe_int_convert(quast.get('num_contigs_0bp')), 'number of contigs >0 bp', 'high')
    
    if not fb:
        fb = ['All tests passed']
    return {'Quality_Module': result, 'Quality_Feedback': '; '.join(fb)}

def process_samples(fastp_files, quast_files, bracken_files=None):
    fastp_map = {extract_sample_id(p): p for p in fastp_files}
    quast_map = {extract_sample_id(p): p for p in quast_files}
    bracken_map = {extract_sample_id(p): p for p in bracken_files} if bracken_files else {}
    all_samples = set(fastp_map.keys()) | set(quast_map.keys()) | set(bracken_map.keys())
    all_data = []
    seen = set()
    
    for sample_id in all_samples:
        if sample_id in seen:
            print(f"Warning: duplicate sample ID {sample_id}")
            continue
        seen.add(sample_id)
        
        entry = {'sample_id': sample_id}
        
        if sample_id in fastp_map:
            try:
                jd = json.load(open(fastp_map[sample_id]))
                entry['fastp'] = parse_fastp_json(jd)
                if not bracken_files:
                    entry['bracken'] = extract_top_taxa_from_bracken(jd)
            except Exception as e:
                print(f"Error parsing FastP for {sample_id}: {e}")
                entry['fastp'] = {}
                if not bracken_files:
                    entry['bracken'] = {}
        else:
            entry['fastp'] = {}
            if not bracken_files:
                entry['bracken'] = {}
        
        if bracken_files:
            if sample_id in bracken_map:
                try:
                    jd = json.load(open(bracken_map[sample_id]))
                    entry['bracken'] = extract_top_taxa_from_bracken(jd)
                except Exception as e:
                    print(f"Error parsing Bracken for {sample_id}: {e}")
                    entry['bracken'] = {}
            else:
                entry['bracken'] = {}
        
        if sample_id in quast_map:
            try:
                entry['quast'] = parse_quast_json(quast_map[sample_id])
            except Exception as e:
                print(f"Error parsing QUAST for {sample_id}: {e}")
                entry['quast'] = {}
        else:
            print(f"Warning: no QUAST file for {sample_id}")
            entry['quast'] = {}
        
        entry['qc'] = calculate_quality_metrics(entry['fastp'], entry['bracken'], entry['quast'])
        all_data.append(entry)
    
    return all_data

def create_aggregated_dataframes(all_data):
    fp_cols = ['sample_id','total_reads','total_bases','q20_rate','q30_rate','gc_content',
               'duplication_rate','passed_filter_reads','low_quality_reads',
               'total_reads_after','total_bases_after','q30_rate_after','gc_content_after','ratio_reads_after']
    br_cols = ['sample_id'] + [f'{p}_{i}' for i in range(1,6) for p in ['taxon','count','proportion']]
    q_keys = set(['sample_id'])
    for e in all_data:
        q_keys.update(e['quast'].keys())
    std = ['N50','Total_length','Largest_contig','num_contigs_500bp','GC_percent']
    qu_cols = ['sample_id'] + [k for k in std if k in q_keys] + sorted([k for k in q_keys if k not in std and k!='sample_id'])
    qc_cols = ['sample_id','Quality_Module','Quality_Feedback',
               'total_reads_after','ratio_reads_after','bracken_main_taxon','proportion_of_main_taxon','reads_gc_content','N50','Total_length','Largest_contig','num_contigs_500bp','num_contigs_0bp']
    
    fp_list, br_list, qu_list, qc_list = [],[],[],[]
    for e in all_data:
        sid = e['sample_id']
        fp = e['fastp']
        br = e['bracken']
        qu = e['quast']
        qc = e['qc']
        
        fp_list.append({c: fp.get(c,np.nan) if c!='sample_id' else sid for c in fp_cols})
        br_list.append({c: br.get(c,np.nan) if c!='sample_id' else sid for c in br_cols})
        qu_list.append({c: qu.get(c,np.nan) if c!='sample_id' else sid for c in qu_cols})
        qc_list.append({
            **{
                'sample_id': sid,
                'Quality_Module': qc['Quality_Module'],
                'Quality_Feedback': qc['Quality_Feedback']
            },
            **{
                'total_reads_after': safe_int_convert(fp.get('total_reads_after')),
                'ratio_reads_after': safe_float_convert(fp.get('ratio_reads_after')),
                'bracken_main_taxon': str(br.get('taxon_1')),
                'proportion_of_main_taxon': safe_float_convert(br.get('proportion_1')),
                'reads_gc_content': safe_float_convert(fp.get('gc_content')),
                'N50': safe_int_convert(qu.get('N50')),
                'Total_length': safe_int_convert(qu.get('Total_length')),
                'Largest_contig': safe_int_convert(qu.get('Largest_contig')),
                'num_contigs_500bp': safe_int_convert(qu.get('num_contigs_500bp')),
                'num_contigs_0bp': safe_int_convert(qu.get('num_contigs_0bp'))
            }
        })
    
    return (
        pd.DataFrame(fp_list, columns=fp_cols),
        pd.DataFrame(qu_list, columns=qu_cols),
        pd.DataFrame(br_list, columns=br_cols),
        pd.DataFrame(qc_list, columns=qc_cols)
    )

def write_enhanced_excel(staramr_xlsx, output_xlsx, fastp_df, quast_df, bracken_df, quality_df):
    try:
        orig = pd.read_excel(staramr_xlsx, sheet_name=None, engine='openpyxl')
    except Exception as e:
        print(f"Error reading {staramr_xlsx}: {e}")
        return
    
    if 'Summary' in orig:
        orig['Summary'] = orig['Summary'].astype({'Quality Module': str,'Quality Module Feedback': str })
        sdf = orig['Summary']
        samp = next((c for c in sdf.columns if 'id' in c.lower()), sdf.columns[0])
        
        fbcol = 'Quality Module'
        if fbcol not in sdf.columns:
            sdf.insert(sdf.columns.get_loc(samp)+1, fbcol, '')
        fb_map = {
            str(r['sample_id']): r['Quality_Module'] for _,r in quality_df.iterrows()
        }
        fb_df = pd.DataFrame(list(fb_map.items()), columns=[samp, fbcol])
        fb_df[samp] = fb_df[samp].astype(str)
        orig['Summary'].update(fb_df, join='left')
        
        fbcol = 'Quality Module Feedback'
        if fbcol not in sdf.columns:
            sdf.insert(sdf.columns.get_loc(samp)+2, fbcol, '')
        fb_map = {
            str(r['sample_id']): r['Quality_Feedback'] for _,r in quality_df.iterrows()
        }
        fb_df = pd.DataFrame(list(fb_map.items()), columns=[samp, fbcol])
        fb_df[samp] = fb_df[samp].astype(str)
        orig['Summary'].update(fb_df, join='left')
        
        fbcol = 'Detected main taxon'
        if fbcol not in sdf.columns:
            sdf.insert(sdf.columns.get_loc(samp)+2, fbcol, '')
        fb_map = {
            str(r['sample_id']): r['bracken_main_taxon'] for _,r in quality_df.iterrows()
        }
        fb_df = pd.DataFrame(list(fb_map.items()), columns=[samp, fbcol])
        fb_df[samp] = fb_df[samp].astype(str)
        orig['Summary'].update(fb_df, join='left')
        
        orig['Summary'].drop(columns=[f'Scheme'], inplace=True)
    else:
        print("Warning: no Summary sheet to update")
    
    for name, df in orig.items():
        if isinstance(df, pd.DataFrame):
            if 'Isolate ID' in df.columns:
                df.sort_values(by='Isolate ID', inplace=True)
                orig[name] = df
    
    if not fastp_df.empty:
        fastp_df.sort_values(by='sample_id', inplace=True)
    if not quast_df.empty:
        quast_df.sort_values(by='sample_id', inplace=True)
    if not bracken_df.empty:
        bracken_df.sort_values(by='sample_id', inplace=True)
    if not quality_df.empty:
        quality_df.sort_values(by='sample_id', inplace=True)
    
    with pd.ExcelWriter(output_xlsx, engine='openpyxl') as w:
        for name, df in orig.items():
            df.to_excel(w, sheet_name=name, index=False)
        if not fastp_df.empty:
            fastp_df.fillna('').to_excel(w, sheet_name='FastP', index=False)
        if not quast_df.empty:
            quast_df.fillna('').to_excel(w, sheet_name='Quast', index=False)
        if not bracken_df.empty:
            bracken_df.fillna('').to_excel(w, sheet_name='Bracken', index=False)
        if not quality_df.empty:
            quality_df.fillna('').to_excel(w, sheet_name='Quality_metrics', index=False)
    
    print(f"Enhanced report saved to {output_xlsx}")

def main():
    p = argparse.ArgumentParser(description="Integrate FastP, Bracken, QUAST into STARAMR Excel report")
    p.add_argument('--staramr', required=True, help="input STARAMR .xlsx")
    p.add_argument('--fastp-jsons', nargs='+', required=True, help="FastP JSON files")
    p.add_argument('--quast-files', nargs='+', required=True, help="QUAST JSON/TSV/CSV files")
    p.add_argument('--bracken-jsons', nargs='+', required=False, help="Bracken JSON files (optional)")
    p.add_argument('--output', required=True, help="output enhanced .xlsx")
    args = p.parse_args()
    
    all_data = process_samples(args.fastp_jsons, args.quast_files, args.bracken_jsons)
    fp_df, qu_df, br_df, qc_df = create_aggregated_dataframes(all_data)
    write_enhanced_excel(args.staramr, args.output, fp_df, qu_df, br_df, qc_df)

if __name__ == '__main__':
    main()
