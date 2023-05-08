import numpy as np
import pandas as pd
import pickle
import os
from pysnptools.snpreader import Bed

from DataGenerator import DataGenerator
from GxEArchitecture import GxEArchitecture
from ld2dependentArchitecture import l2DependentArchitecture
from GxAncestryArchitecture import GxAncestryArchitecture

def DataLoader(params):
    directory = params['directory']
    sample_size = params['sample_size']
    sims = params['sims']
    null_sim = params['null_sim']
    # broad sense heritability
    pve = params['pve']
    # reference data
    use_reference_data = params['use_reference_data']
    sim_win = params['sim_win']
    # parameter for MAF-ES alpha relationship
    exp_alpha = params['exp_alpha']

    snpreader = Bed("/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/data/UK_ALL_REMAINING_1.bed", count_A1=True)
    
    iid = snpreader.iid[:sample_size,1]
    fid = snpreader.iid[:sample_size,0]
    block_size = 5000
    


    if 'ld2' not in params.keys():
        # proportion of linear effects contributing to broad-sense heritability
        rho = params['rho']

        # epistatic sparsity (p_causal)
        sparsity = params['sparsity']
    
    else:
        percentage = params['percentage']

    if 'gxe' in params.keys():
        #n_samples, block_size, pve, rho, psi, win_size, sparsity, phi_sparsity, simtype, exp_alpha = 0.0, null_sim = False
        GxE = GxEArchitecture(sample_size, block_size, pve, rho, params['phi'], params['sim_win'], sparsity, int(params['secondary_sparsity']), params['amplifier'], params['gxe'], exp_alpha)
        y_out, null_alt, linear_pve, epistatic_pve, GxE_pve = GxE.output_data()

        y_out = np.expand_dims(y_out,1)
        null_alt = np.expand_dims(null_alt,1)

        output_prefix = "pheno_"
        output_prefix += "".join([str(np.random.randint(10)) for i in range(10)])        
        output_table = pd.DataFrame(np.array((fid,iid,y_out.flatten())).T,columns=["fid","iid",str(output_prefix)])
        
    
        output_filename = output_prefix + ".txt"

        output_table.to_csv(directory + output_filename, sep="\t", index=False)

        output_obj = {'rho':rho, 'pve':pve, 'win':sim_win, 'sparsity':sparsity, 'secondary_sparsity':params['secondary_sparsity'], 'phi':params['phi'], 'amplifier':params['amplifier'], 'linear_pve':linear_pve, 'epistatic_pve':epistatic_pve, 'gxe_pve':GxE_pve}
        output_filename = output_prefix + ".pkl"

        with open(directory + output_filename, "wb") as f:
            pickle.dump(output_obj, f)

    elif 'ld2' in params.keys():
        print("Generating simulation with ld based architeture")
        #50000, 5000, 0.6, 10000, 10, exp_alpha = 1.0, null_sim = False
        ld2 = l2DependentArchitecture(sample_size, block_size, pve, params['sim_win'], params['percentage'], params['l2sample'],exp_alpha)
        y_out, linear_pve = ld2.output_data()

        y_out = np.expand_dims(y_out,1)
        # null_alt = np.expand_dims(null_alt,1)

        output_prefix = "pheno_"
        output_prefix += "".join([str(np.random.randint(10)) for i in range(10)])        
        output_table = pd.DataFrame(np.array((fid,iid,y_out.flatten())).T,columns=["fid","iid",str(output_prefix)])
        
    
        output_filename = output_prefix + ".txt"

        output_table.to_csv(directory + output_filename, sep="\t", index=False)

        output_obj = {'pve':pve, 'win':sim_win, 'percentage':percentage,  'linear_pve':linear_pve, 'l2sample':params['l2sample']}
        output_filename = output_prefix + ".pkl"

        with open(directory + output_filename, "wb") as f:
            pickle.dump(output_obj, f)
   
    elif 'ancestry' in params.keys():
        print("Generating simulation with GxAncestry architeture")
        #n_samples, block_size, pve, rho, kappa, pc, win_size, sparsity, kappa_sparsity, exp_alpha
        ancestry = GxAncestryArchitecture(sample_size, block_size, pve, rho, params['kappa'], params['pc'], params['sim_win'], sparsity, params['secondary_sparsity'],exp_alpha = 1)
        y_out, null_alt, linear_pve, epistatic_pve, gxancestry_pve = ancestry.output_data()

        y_out = np.expand_dims(y_out,1)
        # null_alt = np.expand_dims(null_alt,1)

        output_prefix = "pheno_"
        output_prefix += "".join([str(np.random.randint(10)) for i in range(10)])        
        output_table = pd.DataFrame(np.array((fid,iid,y_out.flatten())).T,columns=["fid","iid",str(output_prefix)])
        
    
        output_filename = output_prefix + ".txt"

        output_table.to_csv(directory + output_filename, sep="\t", index=False)

        output_obj = {'pve':pve, 'win':sim_win, 'linear_pve':linear_pve, 'pc':params['pc'], 'rho':rho, 'kappa':params['kappa'], 'secondary_sparsity':params['secondary_sparsity'] ,'sparsity':sparsity}
        output_filename = output_prefix + ".pkl"

        with open(directory + output_filename, "wb") as f:
            pickle.dump(output_obj, f)

    else:
        print('DG no GxE')
        DG = DataGenerator(sample_size, block_size, pve, rho, sim_win, sparsity, exp_alpha, null_sim)
        
        y_out, null_alt, linear_pve, epistatic_pve = DG.output_data()
        y_out = np.expand_dims(y_out,1)
        null_alt = np.expand_dims(null_alt,1)

        # output format for ldsc/new meld
        output_prefix = "pheno_"
        output_prefix += "".join([str(np.random.randint(10)) for i in range(10)])
        output_table = pd.DataFrame(np.array((fid,iid,y_out.flatten())).T,columns=["fid","iid",str(output_prefix)])
        

        
        output_filename = output_prefix + ".txt"

        output_table.to_csv(directory + output_filename, sep="\t", index=False)

        output_obj = {'rho':rho, 'pve':pve, 'win':sim_win, 'sparsity':sparsity, 'linear_pve':linear_pve, 'epistatic_pve':epistatic_pve}
        output_filename = output_prefix + ".pkl"

        with open(directory + output_filename, "wb") as f:
            pickle.dump(output_obj, f)
