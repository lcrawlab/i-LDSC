import os
import argparse

from DataDriver import DataLoader

parser = argparse.ArgumentParser()
parser.add_argument("-c","--chi", type=int, default = 50000)
parser.add_argument("-r","--rho", type=float, default = 0.8)
parser.add_argument("-l","--ld", type=int, default = 800)
parser.add_argument("-w","--win", type=int, default = 10000)
parser.add_argument("-p","--pve", type=float, default = 0.6)
parser.add_argument("-s","--sparsity", type=int, default = 10)
parser.add_argument("-d","--ldthreshold", type=float, default = 1.0)
parser.add_argument("-a","--exp-alpha", type=float, default = 0.0)
parser.add_argument("-e","--exppath", type=str, default = "test/")
parser.add_argument("--datagen", dest="datagen", action="store_true")
parser.add_argument("--no-datagen", dest="datagen", action="store_false")
parser.add_argument("--refdata", dest="refdata", action="store_true")
parser.add_argument("--no-refdata", dest="refdata", action="store_false")
#gxe arguments
parser.add_argument("--gxe", dest = "gxe", action = 'store_true')
parser.add_argument("--phi", dest = 'phi', type = float)
parser.add_argument("--ss", dest = 'secondary_sparsity', type = float)
parser.add_argument("--amplifier", dest = 'amplifier', type = float)

#ld score based architecture arguments
parser.add_argument('--ld2based', dest = "ld2", action = 'store_true')
parser.add_argument('--percentage', type=int)
parser.add_argument('--l2sample', dest = 'sample_type')

#ancestry based architecture arguments
parser.add_argument('--ancestry', dest = 'ancestry', action = 'store_true')
parser.add_argument("--kappa", dest = 'kappa', type = float)
parser.add_argument("--pc", dest = 'pc', type = int)

parser.set_defaults(datagen=False)
parser.set_defaults(refdata=False)

args = parser.parse_args()

chi_sample_size = args.chi
ld_sample_size = args.ld
ld_threshold = args.ldthreshold
win_size = args.win
sparsity = args.sparsity
pve = args.pve
datagen = args.datagen
ref_data = args.refdata
exp_path = args.exppath
# relationship between MAF-ES
exp_alpha = args.exp_alpha

if not exp_path.endswith("/"):
    exp_path += "/"
rho = args.rho

if args.gxe:
  phi = args.phi
  secondary_sparsity = args.secondary_sparsity
  amplifier = args.amplifier

if args.ld2:
  percentage = args.percentage
  sample_type = args.sample_type

if args.ancestry:
  kappa = args.kappa
  secondary_sparsity = args.secondary_sparsity
  pc = args.pc
  
"""
ONLY SET PARAMS
"""
params = {
          'sample_size': chi_sample_size, 
          'sims': 1,
          'null_sim': False,
          'ld_sample_size': ld_sample_size,
          'chi_sample_size': chi_sample_size,
          'ld_threshold': ld_threshold,
          'sim_win': win_size,
          'sparsity': sparsity,
          'use_reference_data': ref_data,
          'exp_alpha': exp_alpha
         }

top_directory = "/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/"
sim_directory = top_directory + "sim_cache/" + exp_path
results_directory = top_directory + "sim_results/" + exp_path

params['directory'] = sim_directory
params['sim_directory'] = sim_directory
params['results_directory'] = results_directory
params['experiment'] = exp_path
params['rho'] = rho
params['pve'] = pve

if args.gxe:
  params['phi'] = args.phi
  params['secondary_sparsity'] = args.secondary_sparsity
  params['amplifier'] = args.amplifier
  params['gxe'] = True

if args.ld2:
  params['percentage'] = percentage
  params['l2sample'] = sample_type
  params['ld2'] = True

if args.ancestry:
  params['kappa'] = kappa
  params['secondary_sparsity'] = secondary_sparsity
  params['pc'] = pc
  params['ancestry'] = True


# create directories for this experiment if they do not exist
try:
    if not os.path.exists(sim_directory):
        os.makedirs(sim_directory)
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)
except:
    pass

DataLoader(params)




