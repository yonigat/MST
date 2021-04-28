import numpy as np


def generate_run_name(cfg: dict, ASG: bool=False):
    name = ['M'] if cfg['NUM_PROJECTIONS_IN_SHOT'] > 1 else ['S']
    if ASG:
        name.append('ASG')
    fields = ['NUM_OF_SOURCES', 'NUM_OF_SHOTS_GT', 'NUM_OF_SHOTS_OPTIMIZATION', 'NUM_OF_PHOTONS_FORWARD', 'NUM_OF_PHOTONS_INVERSE']
    for field in fields:
        if type(cfg[field]) is int and cfg[field] > 1000:
            name.append(format_e(cfg[field]))
        else:
            name.append(str(cfg[field]))
    separator = '_'
    full_name = separator.join(name)
    return full_name


def format_e(n):
   a = '%E' % n
   return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]


def parse_run_name(run_name: str):
    if run_name[0] != 'S' and run_name[0] != 'M':
        raise Exception('run name not in desired format')
    run_name = run_name.split(sep='_')
    parsed_name = ['single,'] if run_name[0] == 'S' else ['multiplexed,']
    parsed_name.append('forward photons:')
    parsed_name.append(run_name[4])
    return ' '.join(parsed_name), run_name

def parse_run_name_new(run_name: str):
    if run_name[0] != 'S' and run_name[0] != 'M':
        raise Exception('run name not in desired format')
    parsed_name = {}
    run_name = run_name.split(sep='_')
    parsed_name['multiplexed'] = True if run_name[0] == 'M' else False
    parsed_name['num_sources'] = int(run_name[1])
    parsed_name['num_shots'] = int(run_name[2])
    parsed_name['num_projections_in_shot'] = int(np.ceil(parsed_name['num_sources'] / parsed_name['num_shots']))
    parsed_name['photons_forward'] = float(run_name[4])
    parsed_name['photons_inverse'] = float(run_name[5])
    return parsed_name
