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
      dataset_name = cnv_cmp.split('/')[-1].split('.')[0]
      if dataset_name == 'index':
        continue
      datasets[dataset_name] = cnv_cmp

  out_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'index.json')
  with open(out_path, 'w') as outf:
    print(json.dumps(datasets), file=outf)

main()