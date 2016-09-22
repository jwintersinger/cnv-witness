from __future__ import print_function
import os
import glob
import json
import os
import sys
from collections import defaultdict

def main():
  datasets = {}

  if len(sys.argv) > 1:
    base_dir = os.path.realpath(sys.argv[1])
  else:
    base_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')

  for cnv_cmp in glob.glob(os.path.join(base_dir, '*.json')):
    dataset_name = cnv_cmp.split('/')[-1].split('.')[0]
    if dataset_name in ('index', 'metadata'):
      continue
    with open(cnv_cmp) as cnvf:
      cn_calls = json.load(cnvf)

    datasets[dataset_name] = {
      'cn_calls_path': 'data/%s.json' % dataset_name,
      'methods': cn_calls['methods'],
      'genome_proportions': cn_calls['genome_proportions'],
      'base_cn_state': cn_calls['base_cn_state'],
    }

  out_path = os.path.join(base_dir, 'index.json')
  with open(out_path, 'w') as outf:
    print(json.dumps(datasets), file=outf)

main()
