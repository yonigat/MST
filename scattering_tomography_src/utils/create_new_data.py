from utils.merge_images import merge_images
num_sources_in_shot = [7,8,9] #[2, 3, 4, 5] #[30, 45, 60, 90, 180]#
dir_name = ['knee_'+ str(i) + '_M1E7_NOISY_60' for i in num_sources_in_shot]
single_src_cfg = '../tests/GT_generate_new_cfg.yaml'
for num_source, dir in zip(num_sources_in_shot, dir_name):
    merge_images(single_src_cfg, dir, num_source, non_triv=False)