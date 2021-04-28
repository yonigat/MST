import operator
from utils.lookahead import lookahead
import numpy as np
import itertools

from utils.generate_run_name import parse_run_name_new
from utils.read_TB_logger import read_all_exp_loggers
import matplotlib.pyplot as plt
from utils.general_utils import Extract


def create_iterators():
    color_iter = itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k', 'tab:orange', 'tab:olive', 'tab:brown', 'tab:purple'))
    line_iter = itertools.cycle(('-', '-.', '--', ':', (0, (3, 10, 1, 10))))
    return color_iter, line_iter


exp_dir = '/home/yonatangat/Projects/scattering_tomo/experiments/knee_opt/logs/'

data = read_all_exp_loggers(exp_dir)

eps_ASG = None
delta_ASG = None
for key in list(data):
    if 'ASG' in key:
        eps_ASG = data[key]['epsilon'][1]
        delta_ASG = data[key]['delta'][1]
        data.pop(key, None)

runs = sorted(list(data.keys()))
keys_in_exp = list(data[runs[0]].keys())
keys2plot = ['loss', 'epsilon', 'delta']
num_runs = len(runs)



results = {}


save_num_srcs = True
results['num_srcs'] = []

for key in keys_in_exp:
    if not key in keys2plot:
        continue
    if results.get(key) is None:
        results[key] = []
    for run_ind, run in enumerate(sorted(runs, key=lambda x: (x[0], int(float(x.split(sep='_')[2]))))):
        parsed_run_name = parse_run_name_new(run)
        results[key].append(data[run][key][1][:])
        if save_num_srcs:
            results['num_srcs'].append(parsed_run_name['num_projections_in_shot'])
    save_num_srcs = False

last_loss = Extract(results['loss'])
first_delta = Extract(results['delta'])
last_delta = Extract(results['delta'], ind=-1)
first_epsilon = Extract(results['epsilon'])
last_epsilon = Extract(results['epsilon'], ind=-1)



num_srcs_2_plot = [1, 2, 3, 4, 5]
sorter = np.argsort(results['num_srcs'])
idx_srcs = sorter[np.searchsorted(results['num_srcs'], num_srcs_2_plot, sorter=sorter)]

for key in results.keys():
    colors, lines = create_iterators()
    if key == 'num_srcs':
        continue
    plt.figure()
    # plt.title(f'{key}')
    # plt.ylabel(f'{key}')
    plt.xlabel('iteration')
    if key=='loss':
        pass
        # plt.yscale('log')
    for iter in idx_srcs:
        plt.plot(range(len(results[key][iter])), results[key][iter],
                 linestyle=next(lines), color=next(colors),
                 label=f'{results["num_srcs"][iter]} sources')
    if key=='delta':
        pass
        # plt.ylim(ymax=0)
    else:
        pass
        # plt.ylim(ymin=0)
    # plt.legend()
    plt.show()

# fig = plt.figure()
# plt.title(' relative error in final epsilon compared to single source')
# plt.plot(results['num_srcs'], np.abs(last_epsilon - last_epsilon[-1])/last_delta[-1])
# plt.xlabel('Number of sources multiplexed')
# plt.ylabel('epsilon')
# plt.show()
#
# fig = plt.figure()
# plt.title('relative error in final delta compared to single source')
# plt.plot(results['num_srcs'], np.abs(last_delta - last_delta[-1])/last_delta[-1])
# plt.xlabel('Number of sources multiplexed')
# plt.ylabel('delta')
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111)
# # plt.title('Epsilon vs. number of multiplexed sources')
# ax.plot(results['num_srcs'], last_epsilon, marker='o', label='Scattering Tomography')
# ax.plot(results['num_srcs'], first_epsilon, marker='*', label='Linear Estimation')
# if eps_ASG is not None:
#     ax.scatter(1, eps_ASG, label='with ASG', color='r')
# ax.set_xlabel('Number of sources multiplexed', fontsize=16)
# plt.xticks([1,2,3,4,5], fontsize=16)
# plt.yticks(fontsize=16)
# ax.set_ylabel(r'$\epsilon$')
# plt.legend(fontsize=14)
# plt.ylim(ymin=-0.005)
# plt.show()

## attenuation coefficient insted of density
first_delta = [-0.06985, -0.09422025675420702, -0.10639686102393271, -0.1285091272479237, -0.15115258780890178]
last_delta = [ 0.0032, 0.0019495738910716488, -0.0008523961007311826, -0.0006396413839048884, -0.001669078070629918]
first_epsilon = [0.08376, 0.10205894680340981, 0.11242267019602063, 0.13252064819764767, 0.1542573572365469]
last_epsilon = [0.05387,  0.054692565630227945, 0.05660935865107128, 0.06023317529283833, 0.06600538447029582]


fig = plt.figure(figsize=(14, 14))
ax = fig.add_subplot(111)
# plt.title('Delta vs. number of multiplexed sources')
ax.plot(results['num_srcs'], last_epsilon[::-1], color='b', marker='o', label='Scattering Tomography',markersize=20)
ax.plot(results['num_srcs'], first_epsilon[::-1], color='r', marker='*', label='Linear Estimation', markersize=20)
if delta_ASG is not None:
    ax.scatter(1, delta_ASG, label='with ASG', color='r')
ax.set_xlabel('Number of sources multiplexed', fontsize=32)
plt.xticks([1,2,3,4,5], fontsize=32)
plt.yticks(fontsize=32)
# ax.set_ylabel(r'$\epsilon$', fontsize=28)
plt.legend(fontsize=32)
plt.ylim(ymin=-0.005)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1.0)
plt.show()

fig = plt.figure(figsize=(14, 14))
ax = fig.add_subplot(111)
# plt.title('Delta vs. number of multiplexed sources')
ax.plot(results['num_srcs'], last_delta[::-1], color='b', marker='o', label='Scattering Tomography',markersize=20)
ax.plot(results['num_srcs'], first_delta[::-1],color='r', marker='*', label='Linear Estimation', markersize=20)
if delta_ASG is not None:
    ax.scatter(1, delta_ASG, label='with ASG', color='r')
ax.set_xlabel('Number of sources multiplexed', fontsize=32)
plt.xticks([1,2,3,4,5], fontsize=32)
plt.yticks(fontsize=32)
# ax.set_ylabel(r'$\delta$', fontsize=28)
plt.legend(fontsize=32)
plt.ylim(ymax=0.005)
for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(1.0)
plt.show()


