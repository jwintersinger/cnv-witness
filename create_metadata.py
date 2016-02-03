from __future__ import print_function
import csv
import json
import os
import sys
from collections import defaultdict

def load_tumor_info(tumor_info_fn):
  tumor_info = {}

  with open(tumor_info_fn) as tumorinf:
    reader = csv.DictReader(tumorinf, delimiter='\t')
    for row in reader:
      for sampid in row['tumor_wgs_aliquot_id'].split(','):
        tumor_type = row['dcc_project_code'].split('-')[0]
        tumor_info[sampid] = {
          'tumor_type': tumor_type
        }

  return tumor_info

def load_summaries(summary_fns):
  summaries = defaultdict(dict)

  for summfn in summary_fns:
    method = os.path.basename(summfn).split('.')[0]
    with open(summfn) as summf:
      reader = csv.DictReader(summf, delimiter='\t')
      for row in reader:
        sampid = row['samplename']
        del row['samplename']
        summaries[sampid][method] = row

  return summaries

def main():
  tumor_info_fn = sys.argv[1]
  base_cn_state_fn = sys.argv[2]
  summary_fns = sys.argv[3:]

  tumor_info = load_tumor_info(tumor_info_fn)
  with open(base_cn_state_fn) as basecnf:
    base_cn_state = json.load(basecnf)
  summaries = load_summaries(summary_fns)

  for sampid, base_cn in base_cn_state.items():
    if sampid not in tumor_info:
      continue
    tumor_info[sampid]['base_cn_state'] = base_cn

  for sampid, info in tumor_info.items():
    if sampid not in summaries:
      continue
    tumor_info[sampid]['purity'] = {meth: float(summaries[sampid][meth]['purity']) for meth in summaries[sampid].keys() if summaries[sampid][meth]['purity'] != 'NA'}
    tumor_info[sampid]['ploidy'] = {meth: float(summaries[sampid][meth]['ploidy']) for meth in summaries[sampid].keys() if summaries[sampid][meth]['ploidy'] != 'NA'}

  print(json.dumps(tumor_info))

main()
