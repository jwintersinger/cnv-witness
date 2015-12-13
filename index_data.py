from __future__ import print_function
import os
import glob
import json
import os
from collections import defaultdict

def main():
  datasets = {}

  base_dir = 'data'
  for cnv_cmp in glob.glob(os.path.join(base_dir, '*.json')):
    with open(cnv_cmp) as cnvf:
      cn_calls = json.load(cnvf)

    dataset_name = cnv_cmp.split('/')[-1].split('.')[0]
    if dataset_name in ('index', 'metadata'):
      continue
    datasets[dataset_name] = {
      'cn_calls_path': cnv_cmp,
      'methods': cn_calls['methods'],
      'genome_proportions': cn_calls['genome_proportions']
    }

  out_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'index.json')
  with open(out_path, 'w') as outf:
    print(json.dumps(datasets), file=outf)

main()
