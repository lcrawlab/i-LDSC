import numpy as np
from sklearn.preprocessing import scale
from SNPData import SNPData
import pandas as pd
import sys

class l2DependentArchitecture(object):

    def __init__(self, n_samples, block_size, pve, win_size, perc_cutoff, l2sample, exp_alpha, null_sim = False):
        print("pve {}, top {} percent l2 scores, simtype: {}".format(pve, perc_cutoff, 'ld2 dependent'))
        ## broad sense heritability
        self.pve = pve
        self.perc_cutoff = float(perc_cutoff)/100
        # proportion of linear effects contributing to broad-sense heritability
        chrom_range = np.arange(1,23)
        
        self.l2sample = l2sample
        self.n_samples = n_samples
        self.block_size = block_size
        
        self.snpdata = SNPData(chrom_range, self.n_samples, self.block_size)

        self.null_sim = null_sim
        
        #size of windows that the CASS score is constructed with
        self.window_size = win_size
        
        #Initialize effect size vectors for each type of generative architectures
        self.y_marginal = np.zeros((self.n_samples))
        self.y_err = np.zeros((self.n_samples))


        self.y = np.zeros((self.n_samples))
        # relationship between MAF-ES
        self.exp_alpha = exp_alpha
        
        # for testing how many SNPs make up the cis window
        self.w_shapes = []
        idx = 0

        next_snps = np.zeros((0))
        # counter = 1
        while type(next_snps) == np.ndarray:
            # print(counter)
            if not self.snpdata.curr_chrom_idx >= self.snpdata.n_chromosomes:
                curr_chrom = self.snpdata.chromosomes[self.snpdata.curr_chrom_idx]
                maf = pd.read_csv(f"/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/data/UK_ALL_REMAINING_{curr_chrom}.maf", delim_whitespace = True)

            next_snps, snp_pos, snp_id = self.snpdata.get_next_snps()
            # print(next_snps)
            if type(next_snps) == np.ndarray:
                self.n_snps = next_snps.shape[1]
 
                maf_vals = maf.set_index('SNP').loc[snp_id].MAF.values

                # nan to zero
                next_snps[np.isnan(next_snps)] = 0

                y_marginal, y_err = self.generate_effects(pve, next_snps, snp_id, maf_vals)
                self.y_marginal += y_marginal

                self.y_err += y_err

                idx += 1
            # counter += 1
        self.y_marginal = self.y_marginal*np.sqrt((self.pve)/np.var(self.y_marginal))
        self.y_err=self.y_err*np.sqrt((1-self.pve)/np.var(self.y_err))

        self.y = self.y_marginal + self.y_err
        self.y = scale(self.y)

        self.linear_pve = np.var(self.y_marginal.flatten())/np.var(self.y.flatten())
        
        print("linear pve: ",self.linear_pve)
        print("error pve: ",np.var(self.y_err.flatten())/np.var(self.y.flatten()))

    def generate_effects(self, pve, X, snp_id, maf): 
        #self.pve_per_hub = []
        self.pve = float(pve)
        
        n_snps_draw = {snp_id[i]:i for i in range(self.n_snps)}

        #read in the ld scores and use them to construct the effect size
        l2scores = pd.read_csv('/users/ssmith40/data/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/ldarchitecturescores/whole.genome.l2.ldscore', delim_whitespace = True)
        l2scores = l2scores[['SNP','LDR2_' + str(self.exp_alpha)]]
        
        overlappingdf = l2scores[l2scores['SNP'].isin(snp_id)]
        
        nulls = np.setdiff1d(snp_id,overlappingdf.SNP.tolist())
        
        overlappingdf = overlappingdf.sort_values(by = 'LDR2_' + str(self.exp_alpha))
        
        if self.l2sample == 'high':
            effectfxn_snps = overlappingdf.tail(int(overlappingdf.shape[0]*self.perc_cutoff))
        
        elif self.l2sample =='low':
            effectfxn_snps = overlappingdf.head(int(overlappingdf.shape[0]*self.perc_cutoff))

        noeffectfxn_snps = np.setdiff1d(overlappingdf['SNP'].tolist(),effectfxn_snps['SNP'].tolist())

        nulls = np.concatenate((nulls,noeffectfxn_snps))
        
        effectfxn_snps['genotype_index'] = effectfxn_snps['SNP'].map(n_snps_draw)

        maf = maf[effectfxn_snps['genotype_index'].tolist()]
        self.beta = np.random.normal(0,np.sqrt((2*maf*(1.0-maf))**self.exp_alpha),size = np.shape(maf)[0])
        y_marginal = np.dot(X[:,effectfxn_snps['genotype_index'].tolist()],self.beta)

        y_marginal = y_marginal * np.sqrt(self.pve/(np.var(y_marginal)))


        #sets the error to 1-pve
        # if self.pve > 0:
        #     eps = (1.0-self.pve)*np.var(y_marginal)/self.pve
        # else:
        #     eps = 1
        # y_err = np.random.normal(0,np.sqrt(eps),size=self.n_samples)

        y_err = np.random.normal(0,size=self.n_samples)
        y_err=y_err*np.sqrt((1-self.pve)/np.var(y_err))
    
        y = y_marginal + y_err
        return y_marginal, y_err

    def output_data(self):
        return self.y, self.linear_pve 

#(self, n_samples, block_size, pve, win_size, exp_alpha, null_sim = False)
# x = l2DependentArchitecture(50000, 5000, 0.6, 10000, 10, exp_alpha = 1.0, null_sim = False)
# x