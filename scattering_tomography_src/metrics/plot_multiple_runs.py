import operator

from utils.general_utils import Extract
from utils.lookahead import lookahead

from utils.generate_run_name import parse_run_name, parse_run_name_new, format_e
from utils.read_TB_logger import read_all_exp_loggers
import matplotlib.pyplot as plt

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'tab:brown','tab:orange', 'tab:olive',  'tab:purple']

exp_dir = '/home/yonatangat/Projects/scattering_tomo/experiments/knee_exp_11APR20/logs/'

data = read_all_exp_loggers(exp_dir)

runs = sorted(list(data.keys()))
keys_in_exp = list(data[runs[0]].keys())

results = {'multi': {}, 'single': {}}
results['multi']['num_photons'] = []
results['single']['num_photons'] = []
save_num_photons = True
num_runs = len(runs)

for key in keys_in_exp:
    if results['single'].get(key) is None:
        results['single'][key] = []
        results['multi'][key] = []
    for run_ind, run in enumerate(sorted(runs, key=lambda x: (x[0], int(float(x.split(sep='_')[4]))))):
        parsed_run_name = parse_run_name_new(run)
        # label_run_name = f'{"Multiplexed" if parsed_run_name["multiplexed"] else "Single"}, forward photons: {parsed_run_name["photons_forward"]}'
        # line = '-' if run[0] == 'M' else '--'
        run_type = 'multi' if parsed_run_name['multiplexed'] else 'single'
        results[run_type][key].append(data[run][key][1][:])
        if save_num_photons:
            results[run_type]['num_photons'].append(parsed_run_name['photons_forward'])
    save_num_photons = False

num_photons_plot = [-1, -4, 0]
for key in results['single'].keys():
    if key == 'num_photons':
        continue
    plt.figure()
    plt.title(f'{key}')
    plt.ylabel(f'{key}')
    plt.xlabel('iteration')
    if key=='loss':
        plt.yscale('log')
    for iter in num_photons_plot:
        plt.plot(range(len(results['single'][key][iter])), results['single'][key][iter],
                 linestyle='--', color=colors[iter%len(colors)],
                 label=f'Single, forward photons: {format_e(results["single"]["num_photons"][iter])}')
        plt.plot(range(len(results['multi'][key][iter])), results['multi'][key][iter],
                 linestyle='-.', color=colors[iter%len(colors)],
                 label=f'Multiplexed, forward photons: {format_e(results["multi"]["num_photons"][iter])}')
    plt.legend()
    plt.show()


last_single_loss = Extract(results['single']['loss'])
last_single_delta = Extract(results['single']['delta'])
last_single_epsilon = Extract(results['single']['epsilon'])
last_multi_loss = Extract(results['multi']['loss'])
last_multi_delta = Extract(results['multi']['delta'])
last_multi_epsilon = Extract(results['multi']['epsilon'])

fig = plt.figure()
plt.title('epsilon in multiplexing and single source')
plt.semilogx(results['multi']['num_photons'], last_multi_epsilon, label='Multiplexed')
plt.semilogx(results['single']['num_photons'], last_single_epsilon, '--', label='Single')
plt.xlabel('Number of photons in forward')
plt.ylabel('epsilon')
plt.legend()
plt.show()

fig = plt.figure()
plt.title('delta in multiplexing and single source')
plt.semilogx(results['multi']['num_photons'], last_multi_delta, label='Multiplexed')
plt.semilogx(results['single']['num_photons'], last_single_delta, '--', label='Single')
plt.xlabel('Number of photons in forward')
plt.ylabel('delta')
plt.legend()
plt.show()

fig = plt.figure()
plt.title('final loss in multiplexing and single source')
plt.semilogx(results['multi']['num_photons'], last_multi_loss, label='Multiplexed')
plt.semilogx(results['single']['num_photons'], last_single_loss, '--', label='Single')
plt.xlabel('Number of photons in forward')
plt.ylabel('final loss')
plt.legend()
plt.show()