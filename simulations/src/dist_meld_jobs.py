import os
import sys
import numpy as np
import pickle

exp_path = sys.argv[1] # e.g., fig4
n_jobs = 100
sim_cache_dir = f'/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/sim_cache/{exp_path}/'
results_dir = f'/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/sim_results/{exp_path}/'

filename_prefix = []

for filename in os.listdir(sim_cache_dir):
    if filename.endswith(".txt"):
        filename = filename.rstrip(".txt")
        # for picking up from interrupted results generation
        #if not os.path.isfile(f'{results_dir}{filename}.formatted'):
        #    filename_prefix.append(filename)
        filename_prefix.append(filename)

filename_prefix = np.array(filename_prefix)
filename_prefix = np.array_split(filename_prefix, n_jobs)

job2filename = {}

for idx, filename_set in enumerate(filename_prefix):
    job2filename[idx] = filename_set
# print(job2filename)
pickle.dump(job2filename, open(f'{results_dir}job2filename.{exp_path}.pkl','wb'))

