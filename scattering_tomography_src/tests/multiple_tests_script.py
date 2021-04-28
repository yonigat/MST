import sys
sys.path.append("/home/yonatangat/Projects/scattering_tomo/scattering_tomography_src/")
import os
import numpy as np
from utils.delete_files import clear_dir
import yaml
import subprocess

tests = []
params = {}

params['NUM_PROJECTIONS_IN_SHOT'] = [4] #[1, 2, 4, 5]
params['NUM_OF_SOURCES'] = [180 for _ in params['NUM_PROJECTIONS_IN_SHOT']]
params['NUM_OF_SHOTS_GT'] = [int(x/y) for x, y in zip(params['NUM_OF_SOURCES'], params['NUM_PROJECTIONS_IN_SHOT'])]
# params['NUM_OF_PHOTONS_FORWARD'] = [int(1e6/x) for x in params['NUM_PROJECTIONS_IN_SHOT']]
params['NUM_OF_PHOTONS_FORWARD_ORIG'] = [int(5e7/x) for x in params['NUM_PROJECTIONS_IN_SHOT']]
params['NUM_OF_PHOTONS_INVERSE'] = [int(1e6/x) for x in params['NUM_PROJECTIONS_IN_SHOT']]
params['OPTIMIZATION_ITER'] = [50] #[50, 75, 100, 100]
params['GT_DATA_DIR'] = ['knee_'+str(i)+'_M1E8_NEW' for i in params['NUM_PROJECTIONS_IN_SHOT']]
params['DATA_DIR'] = params['GT_DATA_DIR']

# params['NUM_OF_PHOTONS_FORWARD'] = [int(num) for num in [5e2, 1e3, 3e3, 5e3, 1e4, 1e5, 5e5, 1e6, 3e6]]
# params['NUM_PROJECTIONS_IN_SHOT'] = [int(x/y) for x, y in zip(params['NUM_OF_SOURCES'], params['NUM_OF_SHOTS_GT'])]
#
# params['DATA_DIR'] = ['knee_S5E2', 'knee_S1E3', 'knee_S3E3', 'knee_S5E3', 'knee_S1E4', 'knee_S1E5', 'knee_S5E5',
# 'knee_S1E6', 'knee_S3E6']
# params['DATA_DIR'] = ['knee_M5E2', 'knee_M1E3', 'knee_M3E3', 'knee_M5E3', 'knee_M1E4',
# 'knee_M1E5', 'knee_M5E5', 'knee_M1E6', 'knee_M3E6']


exp_dir = 'density_adjust'

basic_yaml_file = 'GT_generate_new_cfg.yaml'
updated_yaml_file = 'full_optimization_cfg.yaml'
script_to_run = 'full_optimization_script.py'

basic_yaml_file = os.path.join(os.path.abspath(''), basic_yaml_file)
updated_yaml_file = os.path.join(os.path.abspath(''), updated_yaml_file)
script_to_run = os.path.join(os.path.abspath(''), script_to_run)

for _ in range(len(params[next(iter(params))])):
    tests.append({})

for key, vals in params.items():
    for ind, val in enumerate(vals):
        tests[ind][key] = val

with open(basic_yaml_file) as f:
    basic_cfg = yaml.safe_load(f)

for test in tests:
    curr_cfg = basic_cfg
    curr_cfg['EXP_DIR'] = exp_dir
    for key, val in test.items():
        curr_cfg[key] = val
    with open(updated_yaml_file, 'w') as f:
        yaml.safe_dump(curr_cfg, f)
    p = subprocess.Popen(['python', script_to_run])
    p.wait()

